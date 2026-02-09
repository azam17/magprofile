/*
 * gui_analysis.h — Background analysis pipeline for HalalSeq GUI.
 *
 * Runs index loading, FASTQ reading, classification, EM, and report
 * generation on a background thread.  The main (GUI) thread polls the
 * state/progress fields to update the UI.
 */
#ifndef HALALSEQ_GUI_ANALYSIS_H
#define HALALSEQ_GUI_ANALYSIS_H

#include "report.h"
#include <SDL.h>

typedef enum {
    ANALYSIS_IDLE,
    ANALYSIS_LOADING_INDEX,
    ANALYSIS_READING_FASTQ,
    ANALYSIS_CLASSIFYING,
    ANALYSIS_RUNNING_EM,
    ANALYSIS_GENERATING_REPORT,
    ANALYSIS_DONE,
    ANALYSIS_ERROR
} analysis_state_t;

typedef struct {
    /* ---- inputs (set by GUI thread before analysis_start) ---- */
    char fastq_path[1024];
    char index_path[1024];

    /* ---- observable state (read by GUI thread, written by worker) ---- */
    volatile analysis_state_t state;
    volatile int progress_reads;      /* reads processed so far            */
    volatile int progress_total;      /* total reads (0 until known)       */
    char error_msg[256];

    /* ---- output (valid once state == ANALYSIS_DONE) ---- */
    halal_report_t *report;           /* NULL until done                   */

    /* ---- internal ---- */
    SDL_Thread *thread;
} analysis_context_t;

/* Initialise to all zeros / ANALYSIS_IDLE. */
void analysis_init(analysis_context_t *ctx);

/* Launch the background worker.  ctx->fastq_path and ctx->index_path
   must already be set.  Returns 0 on success, -1 on error. */
int  analysis_start(analysis_context_t *ctx);

/* Free any resources held by a previous run (report, thread handle).
   Safe to call multiple times or on a freshly-init'd context. */
void analysis_cleanup(analysis_context_t *ctx);

/* State label strings for the UI status line. */
const char *analysis_state_label(analysis_state_t s);

#endif /* HALALSEQ_GUI_ANALYSIS_H */
