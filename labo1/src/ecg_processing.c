#include <ecg_processing.h>
#include <output_structs.h>
#include <stddef.h>
#include <ecg_utils.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


ECG_Context *ecg_create(const ECG_Params *params) {
    if(!params)
        return NULL;
    ECG_Context* ctx = malloc(sizeof(ECG_Context));
    if (ctx == NULL) {
        return NULL;
    }
    ctx->current_data = malloc(sizeof(double) * MAX_SAMPLES);
    if (ctx->current_data == NULL) {
        free(ctx);
        return NULL;
    }
    ctx->processed_data = malloc(sizeof(double) * MAX_SAMPLES);
    if (ctx->processed_data == NULL) {
        free(ctx->current_data);
        free(ctx);
        return NULL;
    }

    ctx->original_signal = malloc(sizeof(double) * MAX_SAMPLES);
    if (ctx->original_signal == NULL) {
        free(ctx->processed_data);
        free(ctx->current_data);
        free(ctx);
        return NULL;
    }

    ctx->params = params;
    ctx->current_state = ECG_BEFORE_LOWPASS;
    ctx->data_size = MAX_SAMPLES;

    // Initialiser les temps à 0
    memset(&ctx->timing, 0, sizeof(ECG_Timing));

    return ctx;
}


void ecg_destroy(ECG_Context *ctx) {
    ctx->params = NULL;
    free(ctx->current_data);
    free(ctx->processed_data);
    free(ctx->original_signal);
    free(ctx);
}


ECG_Status ecg_analyze(
    ECG_Context *ctx,
    const double *signal,
    size_t n_samples,
    int lead_idx,
    ECG_Peaks *peaks,
    ECG_Intervals *intervals) {

    clock_t start_total = clock();
/**
    * étapes a faire :
    1. Filtrage passe-bande (environ 5–15 Hz),
    2. Calcul de la dérivée du signal,
    3. Mise au carré (accentuation des pentes),
    4. Intégration sur fenêtre glissante,
    5. Détection de seuil adaptatif.
 **/
    if (!ctx || !signal || !peaks) return ECG_ERR_NULL;
    if (lead_idx < 0 || lead_idx >= LEADS) return ECG_ERR_PARAM;
    if (n_samples == 0 || n_samples > MAX_SAMPLES) return ECG_ERR_PARAM;

    memset(peaks, 0, sizeof(ECG_Peaks));
    if (intervals) {
        memset(intervals, 0, sizeof(ECG_Intervals));
    }

    ctx->data_size = n_samples;
    memcpy(ctx->current_data, signal, sizeof(double) * n_samples);
    memcpy(ctx->original_signal, signal, sizeof(double) * n_samples);

// remove le offset bizzare
    ecg_remove_dc(ctx->current_data, n_samples);
    ecg_remove_dc(ctx->original_signal, n_samples);

// check le gain et l'appliquer si nécessaire
    if (ctx->params->gain != 0.0 && ctx->params->gain != 1.0) {
        ecg_apply_gain(ctx->current_data, n_samples, ctx->params->gain);
        ecg_apply_gain(ctx->original_signal, n_samples, ctx->params->gain);
    }
// 1. Filtrage passe-bande (environ 5–15 Hz)

    // low pass
    clock_t low_pass_start = clock();

    if (!ctx || !ctx->current_data || !ctx->processed_data) {
        ctx->current_state = ECG_BEFORE_LOWPASS;
        return ECG_ERR_FAIL;
    }

    ecg_moving_average(ctx->current_data, ctx->processed_data, ctx->data_size, ctx->params->lowpass_window ? ctx->params->lowpass_window : 1);

    swap(ctx->current_data, ctx->processed_data);
    ctx->current_state = ECG_BEFORE_HIGPASS;

    clock_t low_pass_end = clock();
    ctx->timing.lowpass_time = (double)(low_pass_end - low_pass_start) / CLOCKS_PER_SEC;


    // high pass
    clock_t high_pass_start = clock();


    if (!ctx || !ctx->current_data || !ctx->processed_data) {
        ctx->current_state = ECG_BEFORE_HIGPASS;
        return ECG_ERR_FAIL;
    }
    ecg_highpass_ma(ctx->current_data, ctx->processed_data, ctx->data_size, ctx->params->highpass_window ? ctx->params->highpass_window : 1);

    swap(ctx->current_data, ctx->processed_data);
    ctx->current_state = ECG_BEFORE_DERIVATIVE;

    clock_t high_pass_end = clock();
    ctx->timing.highpass_time = (double)(high_pass_end - high_pass_start) / CLOCKS_PER_SEC;

// 2. Calcul de la dérivée du signal,
    clock_t derivative_start = clock();

    if (!ctx || !ctx->current_data || !ctx->processed_data) {
        ctx->current_state = ECG_BEFORE_DERIVATIVE;
        return ECG_ERR_FAIL;
    }
    ecg_derivative_1(ctx->current_data, ctx->processed_data, ctx->data_size);

    swap(ctx->current_data, ctx->processed_data);
    ctx->current_state = ECG_BEFORE_SQUARE;

    clock_t derivative_end = clock();
    ctx->timing.derivative_time = (double)(derivative_end - derivative_start) / CLOCKS_PER_SEC;

// 3. Mise au carré (accentuation des pentes),
    clock_t square_start = clock();

    if (!ctx || !ctx->current_data || !ctx->processed_data) {
        ctx->current_state = ECG_BEFORE_SQUARE;
        return ECG_ERR_FAIL;
    }
    ecg_square(ctx->current_data, ctx->processed_data, ctx->data_size);

    swap(ctx->current_data, ctx->processed_data);
    ctx->current_state = ECG_BEFORE_MWI;

    clock_t square_end = clock();
    ctx->timing.square_time = (double)(square_end - square_start) / CLOCKS_PER_SEC;


// 4. Intégration sur fenêtre glissante,he number of samples to average is chosen in order to average on windows of 150 ms.
    clock_t mwi_start = clock();

    if (!ctx || !ctx->current_data || !ctx->processed_data) {
        ctx->current_state = ECG_BEFORE_MWI;
        return ECG_ERR_FAIL;
    }

    ecg_mwi(ctx->current_data, ctx->processed_data, ctx->data_size, ctx->params->mwi_window? ctx->params->mwi_window : 1);

    swap(ctx->current_data, ctx->processed_data);
    ctx->current_state = ECG_BEFORE_DECISION;

    clock_t mwi_end = clock();
    ctx->timing.mwi_time = (double)(mwi_end - mwi_start) / CLOCKS_PER_SEC;


// 5 détecter les pick

    clock_t peaks_start = clock();

    if (!ctx || !ctx->current_data || !peaks) {
        ctx->current_state = ECG_BEFORE_DECISION;
        return ECG_ERR_FAIL;
    }

    peaks->R_count = 0;

    double max_val = 0.0;
    for (size_t i = 0; i < ctx->data_size; i++)
        if (ctx->current_data[i] > max_val)
            max_val = ctx->current_data[i];

    double threshold = 0.6 * max_val;

    size_t refractory = ctx->params->sampling_rate_hz / 5;
    size_t next = 0;

    for (size_t i = 0; i < ctx->data_size; i++) {

        /* condition maximum local */
        if (ctx->current_data[i] > ctx->current_data[i - 1] &&
            ctx->current_data[i] > ctx->current_data[i + 1] &&
            ctx->current_data[i] > threshold &&
            i > next)
        {
            // Look in window around the detected peak in the original signal to find the true R peak (since we are on the MWI signal which is delayed and smoothed)
            size_t search_start = (i > refractory / 2) ? (i - refractory / 2) : 0;
            size_t search_end = (i + refractory / 2 < ctx->data_size) ? (i + refractory / 2) : ctx->data_size - 1;
            size_t true_peak = i;
            for (size_t j = search_start; j <= search_end; j++) {
                if (ctx->original_signal[j] > ctx->original_signal[true_peak]) {
                    true_peak = j;
                }
            }
            if (peaks->R_count < MAX_BEATS) {
                peaks->R[peaks->R_count++] = true_peak;
            }
            next = i + refractory;
        }
    }

    ctx->current_state = ECG_FINISHED;

    clock_t peaks_end = clock();
    ctx->timing.peaks_time = (double)(peaks_end - peaks_start) / CLOCKS_PER_SEC;


    // intervalle RR
    if (intervals && peaks->R_count > 1) {
        intervals->count = peaks->R_count - 1;
        for (int i = 0; i < intervals->count; i++) {
            double rr_samples = (double)(peaks->R[i+1] - peaks->R[i]);
            intervals->RR[i] = rr_samples / (double)ctx->params->sampling_rate_hz;
        }
    }

    clock_t end_total = clock();
    ctx->timing.total_time = (double)(end_total - start_total) / CLOCKS_PER_SEC;

    return ECG_OK;
}
