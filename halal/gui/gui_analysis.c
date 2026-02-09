/*
 * gui_analysis.c — Background analysis pipeline for HalalSeq GUI.
 *
 * Mirrors the logic in src/main.c cmd_run(), but runs on an SDL thread
 * and reports incremental state/progress to the GUI.
 */
#include "gui_analysis.h"
#include "utils.h"
#include "refdb.h"
#include "index.h"
#include "fastq.h"
#include "classify.h"
#include "em.h"
#include "report.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* State labels                                                        */
/* ------------------------------------------------------------------ */
static const char *state_labels[] = {
    [ANALYSIS_IDLE]              = "Ready",
    [ANALYSIS_LOADING_INDEX]     = "Loading index...",
    [ANALYSIS_READING_FASTQ]     = "Reading FASTQ...",
    [ANALYSIS_CLASSIFYING]       = "Classifying reads...",
    [ANALYSIS_RUNNING_EM]        = "Running EM...",
    [ANALYSIS_GENERATING_REPORT] = "Generating report...",
    [ANALYSIS_DONE]              = "Done",
    [ANALYSIS_ERROR]             = "Error",
};

const char *analysis_state_label(analysis_state_t s) {
    if (s < 0 || s > ANALYSIS_ERROR) return "Unknown";
    return state_labels[s];
}

/* ------------------------------------------------------------------ */
/* Worker thread                                                       */
/* ------------------------------------------------------------------ */
static int analysis_worker(void *data) {
    analysis_context_t *ctx = (analysis_context_t *)data;

    /* 1. Load index ------------------------------------------------- */
    ctx->state = ANALYSIS_LOADING_INDEX;
    halal_index_t *idx = index_load(ctx->index_path);
    if (!idx) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to load index: %s", ctx->index_path);
        ctx->state = ANALYSIS_ERROR;
        return 1;
    }

    /* 2. Read FASTQ ------------------------------------------------- */
    ctx->state = ANALYSIS_READING_FASTQ;
    char **seqs = NULL;
    char **names = NULL;
    int  *lens = NULL;
    int   n_reads = 0;

    if (hs_fasta_read_all(ctx->fastq_path, &seqs, &names, &lens, &n_reads) < 0) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to read: %s", ctx->fastq_path);
        index_destroy(idx);
        ctx->state = ANALYSIS_ERROR;
        return 1;
    }
    ctx->progress_total = n_reads;
    ctx->progress_reads = n_reads;  /* reads are already loaded */

    /* 3. Classify --------------------------------------------------- */
    ctx->state = ANALYSIS_CLASSIFYING;
    classify_opts_t copts = classify_opts_default();
    read_result_t *results = classify_reads(idx, (const char **)seqs, lens,
                                             n_reads, &copts);
    if (!results) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Classification failed");
        hs_fasta_free_all(seqs, names, lens, n_reads);
        index_destroy(idx);
        ctx->state = ANALYSIS_ERROR;
        return 1;
    }

    /* 4. EM --------------------------------------------------------- */
    ctx->state = ANALYSIS_RUNNING_EM;
    int n_em_reads = 0;
    em_read_t *em_reads = em_reads_from_classify(results, n_reads, &n_em_reads);

    em_config_t ecfg = em_config_default();

    /* Mito copy number priors */
    double *mito_cn = (double *)hs_malloc(
        (size_t)idx->db->n_species * sizeof(double));
    for (int s = 0; s < idx->db->n_species; s++)
        mito_cn[s] = idx->db->species[s].mito_copy_number;
    ecfg.mito_copy_numbers = mito_cn;

    /* Amplicon lengths */
    int *amp_lens = (int *)hs_calloc(
        (size_t)(idx->db->n_species * idx->db->n_markers), sizeof(int));
    for (int i = 0; i < idx->db->n_marker_refs; i++) {
        int si = idx->db->markers[i].species_idx;
        int mi = idx->db->markers[i].marker_idx;
        amp_lens[si * idx->db->n_markers + mi] =
            idx->db->markers[i].amplicon_length;
    }

    em_result_t *em = NULL;
    if (n_em_reads > 0) {
        em = em_fit(em_reads, n_em_reads,
                    idx->db->n_species, idx->db->n_markers,
                    amp_lens, &ecfg);
    }

    /* 5. Report ----------------------------------------------------- */
    ctx->state = ANALYSIS_GENERATING_REPORT;
    halal_report_t *report = NULL;
    if (em) {
        report = report_generate(em, idx->db, results, n_reads, 0.001);
    } else {
        /* No classified reads — generate a minimal report */
        report = (halal_report_t *)hs_calloc(1, sizeof(halal_report_t));
        snprintf(report->sample_id, sizeof(report->sample_id), "%s",
                 ctx->fastq_path);
        report->verdict = INCONCLUSIVE;
        report->total_reads = n_reads;
        report->classified_reads = 0;
    }

    /* Cleanup intermediaries */
    if (em) em_result_destroy(em);
    em_reads_free(em_reads, n_em_reads);
    classify_results_free(results, n_reads);
    free(amp_lens);
    free(mito_cn);
    hs_fasta_free_all(seqs, names, lens, n_reads);
    index_destroy(idx);

    /* Publish result */
    ctx->report = report;
    ctx->state = ANALYSIS_DONE;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */
void analysis_init(analysis_context_t *ctx) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->state = ANALYSIS_IDLE;
}

int analysis_start(analysis_context_t *ctx) {
    /* Clean up any previous run */
    analysis_cleanup(ctx);

    ctx->state        = ANALYSIS_IDLE;
    ctx->report       = NULL;
    ctx->progress_reads = 0;
    ctx->progress_total = 0;
    ctx->error_msg[0]   = '\0';

    ctx->thread = SDL_CreateThread(analysis_worker, "analysis", ctx);
    if (!ctx->thread) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to create thread: %s", SDL_GetError());
        ctx->state = ANALYSIS_ERROR;
        return -1;
    }
    return 0;
}

void analysis_cleanup(analysis_context_t *ctx) {
    if (ctx->thread) {
        int status;
        SDL_WaitThread(ctx->thread, &status);
        ctx->thread = NULL;
    }
    if (ctx->report) {
        report_destroy(ctx->report);
        ctx->report = NULL;
    }
}
