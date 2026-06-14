#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>


#include "csv_reader.h"
#include "ecg_processing.h"
#include "json_writer.h"
#include "output_structs.h"

#define PACKET_SIZE 1000
#define OVERLAP_SIZE 250

static void append_r_peak(ECG_Peaks *global_peaks, int sample_index)
{
    if (global_peaks->R_count < MAX_BEATS) {
        global_peaks->R[global_peaks->R_count++] = sample_index;
    }
}

static void append_rr_interval(ECG_Intervals *global_intervals, double rr_value)
{
    if (global_intervals->count < MAX_BEATS) {
        global_intervals->RR[global_intervals->count++] = rr_value;
    }
}

// Structure pour stocker les résultats d'un paquet
typedef struct {
    ECG_Peaks   peaks;
    ECG_Intervals intervals;
    ECG_Timing  timing;
    size_t      packet_start;
    size_t      accept_until;
    int         valid;
} PacketResult;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_csv> <output_json>\n", argv[0]);
        return 1;
    }

    if (read_csv(argv[1]) != 0) {
        fprintf(stderr, "Erreur lecture CSV.\n");
        return 2;
    }

    ECG_Peaks peaks;
    ECG_Intervals intervals;
    memset(&peaks, 0, sizeof(peaks));
    memset(&intervals, 0, sizeof(intervals));

    ECG_Params params;
    memset(&params, 0, sizeof(params));
    params.sampling_rate_hz = SAMPLING_RATE;
    params.leads            = LEADS;
    params.gain             = 1.0; // <-- Ajuster le gain si nécessaire
    params.r_threshold_hint = 0.0; // <-- Optionnel, peut être 0.0, et peut-être adaptatif au long du code.

    params.dc_offset        = 13.0;

    // Paramètres de l'algorithme de Pan-Tompkins
    // Fenêtre du filtre passe-bas : environ 30-35 ms à 500 Hz = 15-17 échantillons
    params.lowpass_window   = 15;

    // Fenêtre du filtre passe-haut : environ 200-400 ms à 500 Hz = 100-200 échantillons
    params.highpass_window  = 125;

    // Fenêtre d'intégration (MWI) : environ 150 ms à 500 Hz = 75 échantillons
    params.mwi_window       = 75;

    params.r_threshold_hint = 0.0; // Seuil adaptatif (calculé automatiquement)


    ECG_Context *ctx = ecg_create(&params);
    if (!ctx) {
        fprintf(stderr, "Erreur: ecg_create() a échoué.\n");
        return 4;
    }

    int lead_index = 1; // Analyser la LEAD II (index 1)
    if (lead_index < 0 || lead_index >= LEADS) {
        fprintf(stderr, "Erreur: lead_index invalide.\n");
        ecg_destroy(ctx);
        return 5;
    }

    clock_t start = clock();

    const size_t packet_size = PACKET_SIZE;
    const size_t overlap_size = OVERLAP_SIZE;
    const size_t stride = packet_size - overlap_size;

    if (stride == 0) {
        fprintf(stderr, "Erreur: chevauchement invalide.\n");
        ecg_destroy(ctx);
        return 7;
    }

    // Calculer le nombre total de paquets
    size_t num_packets = 0;
    for (size_t s = 0; s < (size_t)sample_count; s += stride) {
        num_packets++;
        if (s + packet_size >= (size_t)sample_count) break;
    }

    // Allouer les résultats et les contextes par thread
    PacketResult *results = calloc(num_packets, sizeof(PacketResult));
    if (!results) { /* erreur */ return 8; }

    ECG_Timing total_timing;
    memset(&total_timing, 0, sizeof(total_timing));

    #pragma omp parallel
    {
        // Chaque thread a son propre contexte ECG
        ECG_Context *local_ctx = ecg_create(&params);

        #pragma omp for schedule(dynamic, 1)
        for (size_t p = 0; p < num_packets; p++) {
            size_t packet_start = p * stride;
            size_t packet_end   = packet_start + packet_size;
            if (packet_end > (size_t)sample_count)
                packet_end = (size_t)sample_count;

            size_t packet_samples = packet_end - packet_start;

            PacketResult *res = &results[p];
            res->packet_start = packet_start;
            res->accept_until = packet_start + stride;
            if (res->accept_until > (size_t)sample_count)
                res->accept_until = (size_t)sample_count;

            memset(&res->peaks,     0, sizeof(res->peaks));
            memset(&res->intervals, 0, sizeof(res->intervals));

            ECG_Status st = ecg_analyze(
                local_ctx,
                &ecg_data[lead_index][packet_start],
                packet_samples,
                lead_index,
                &res->peaks,
                &res->intervals
            );

            if (st == ECG_OK) {
                res->valid = 1;
                res->timing = local_ctx->timing;
            }
        }

        ecg_destroy(local_ctx);

        // Accumulation thread-safe des timings
        #pragma omp critical
        {
            for (size_t p = 0; p < num_packets; p++) {
                if (!results[p].valid) continue;
                total_timing.lowpass_time    += results[p].timing.lowpass_time;
                total_timing.highpass_time   += results[p].timing.highpass_time;
                total_timing.derivative_time += results[p].timing.derivative_time;
                total_timing.square_time     += results[p].timing.square_time;
                total_timing.mwi_time        += results[p].timing.mwi_time;
                total_timing.peaks_time      += results[p].timing.peaks_time;
                total_timing.total_time      += results[p].timing.total_time;
            }
        }
    }


    // Réassemblage séquentiel dans l'ordre des paquets
    for (size_t p = 0; p < num_packets; p++) {
        PacketResult *res = &results[p];
        if (!res->valid) continue;

        for (int i = 0; i < res->peaks.R_count; i++) {
            size_t global_index = res->packet_start + (size_t)res->peaks.R[i];
            if (global_index < res->accept_until) {
                append_r_peak(&peaks, (int)global_index);
            }
        }
    }

    clock_t end = clock();
    double elapsed_sec = (double)(end - start) / CLOCKS_PER_SEC;
    free(results);

    printf("Analyse ECG: %.6f s\n", elapsed_sec);

    if (peaks.R_count > 1) {
        intervals.count = 0;
        for (int i = 0; i < peaks.R_count - 1; i++) {
            double rr_samples = (double)(peaks.R[i + 1] - peaks.R[i]);
            append_rr_interval(&intervals, rr_samples / (double)ctx->params->sampling_rate_hz);
        }
    }

    printf("%d pics R détectés.\n", peaks.R_count);
    printf("\n=== Performance Analysis ===\n");
    printf("low pass   : %.6f s\n", total_timing.lowpass_time);
    printf("high pass  : %.6f s\n", total_timing.highpass_time);
    printf("Derivative : %.6f s\n", total_timing.derivative_time);
    printf("squarre    : %.6f s\n", total_timing.square_time);
    printf("mwi        : %.6f s\n", total_timing.mwi_time);
    printf("peaks      : %.6f s\n", total_timing.peaks_time);
    printf("---------------------------------------------\n");
    printf("Total      : %.6f s\n", total_timing.total_time);
    ecg_destroy(ctx);

    if (write_json(argv[2], &peaks, &intervals) != 0) {
        fprintf(stderr, "Erreur écriture JSON.\n");
        return 3;
    }

    printf("Analyse terminée. Résultats sauvegardés dans %s\n", argv[2]);
    return 0;
}
