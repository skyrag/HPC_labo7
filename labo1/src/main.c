#include <stdio.h>
#include <string.h>
#include <time.h>

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

    ECG_Timing total_timing;
    memset(&total_timing, 0, sizeof(total_timing));

    const size_t packet_size = PACKET_SIZE;
    const size_t overlap_size = OVERLAP_SIZE;
    const size_t stride = packet_size - overlap_size;

    if (stride == 0) {
        fprintf(stderr, "Erreur: chevauchement invalide.\n");
        ecg_destroy(ctx);
        return 7;
    }

    for (size_t packet_start = 0; packet_start < (size_t)sample_count; packet_start += stride) {
        size_t packet_end = packet_start + packet_size;
        if (packet_end > (size_t)sample_count) {
            packet_end = (size_t)sample_count;
        }

        size_t packet_samples = packet_end - packet_start;
        if (packet_samples == 0) {
            break;
        }

        ECG_Peaks packet_peaks;
        ECG_Intervals packet_intervals;
        memset(&packet_peaks, 0, sizeof(packet_peaks));
        memset(&packet_intervals, 0, sizeof(packet_intervals));

        ECG_Status st = ecg_analyze(
            ctx,
            &ecg_data[lead_index][packet_start],
            packet_samples,
            lead_index,
            &packet_peaks,
            &packet_intervals
        );

        total_timing.lowpass_time += ctx->timing.lowpass_time;
        total_timing.highpass_time += ctx->timing.highpass_time;
        total_timing.derivative_time += ctx->timing.derivative_time;
        total_timing.square_time += ctx->timing.square_time;
        total_timing.mwi_time += ctx->timing.mwi_time;
        total_timing.peaks_time += ctx->timing.peaks_time;
        total_timing.total_time += ctx->timing.total_time;

        if (st != ECG_OK) {
            fprintf(stderr, "Erreur: ecg_analyze() a retourné %d sur le paquet commençant à %zu.\n", (int)st, packet_start);
            ecg_destroy(ctx);
            return 6;
        }

        size_t accept_until = packet_start + stride;
        if (accept_until > (size_t)sample_count) {
            accept_until = (size_t)sample_count;
        }

        for (int i = 0; i < packet_peaks.R_count; i++) {
            size_t global_index = packet_start + (size_t)packet_peaks.R[i];
            if (global_index < accept_until) {
                append_r_peak(&peaks, (int)global_index);
            }
        }

        if (packet_start + packet_size >= (size_t)sample_count) {
            break;
        }
    }

    clock_t end = clock();
    double elapsed_sec = (double)(end - start) / CLOCKS_PER_SEC;

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
