/*
 * main_gui.c — HalalSeq desktop GUI application.
 *
 * SDL2 + Nuklear immediate-mode GUI.  Single window (960x640), left
 * panel for file selection / progress, right panel for results.
 *
 * Build:  see Makefile "gui" target or CMakeLists.txt.
 */

/* ================================================================== */
/* Nuklear implementation (must come first, once per translation unit) */
/* ================================================================== */
#include "nuklear_setup.h"

/* ================================================================== */
/* Project headers                                                     */
/* ================================================================== */
#include "gui_analysis.h"
#include "tinyfiledialogs.h"
#include "refdb.h"

#include <stdio.h>
#include <string.h>

/* ================================================================== */
/* Constants                                                           */
/* ================================================================== */
#define WINDOW_W  960
#define WINDOW_H  640
#define LEFT_W    340
#define VERSION   "0.1.0"

/* Row heights (logical points) */
#define ROW_TITLE  30
#define ROW_LABEL  26
#define ROW_BTN    38
#define ROW_SMALL  22
#define ROW_TABLE  24
#define ROW_SPACE  8
#define ROW_BAR    14

/* ================================================================== */
/* Colour helpers                                                      */
/* ================================================================== */
static struct nk_color col_pass       = {34, 139, 34, 255};   /* forest green */
static struct nk_color col_fail       = {220, 20, 60, 255};   /* crimson      */
static struct nk_color col_inconc     = {218, 165, 32, 255};  /* goldenrod    */
static struct nk_color col_halal      = {34, 139, 34, 255};
static struct nk_color col_haram      = {220, 20, 60, 255};
static struct nk_color col_mashbooh   = {218, 165, 32, 255};
static struct nk_color col_unknown    = {160, 160, 160, 255};

static struct nk_color status_color(halal_status_t s) {
    switch (s) {
        case HALAL:            return col_halal;
        case HARAM:            return col_haram;
        case MASHBOOH:         return col_mashbooh;
        default:               return col_unknown;
    }
}

static struct nk_color verdict_color(verdict_t v) {
    switch (v) {
        case PASS:         return col_pass;
        case FAIL:         return col_fail;
        default:           return col_inconc;
    }
}

/* ================================================================== */
/* Friendly species name lookup                                        */
/* ================================================================== */
static const char *friendly_species_name(const char *species_id) {
    /* Map Latin names to common food names */
    struct { const char *latin; const char *common; } names[] = {
        { "Bos_taurus",       "Beef (Cow)"           },
        { "Sus_scrofa",       "Pork (Pig)"           },
        { "Ovis_aries",       "Lamb (Sheep)"         },
        { "Gallus_gallus",    "Chicken"              },
        { "Capra_hircus",     "Goat"                 },
        { "Equus_caballus",   "Horse"                },
        { "Bubalus_bubalis",  "Buffalo"              },
        { "Anas_platyrhynchos","Duck"                },
        { "Cervus_elaphus",   "Deer (Venison)"       },
        { "Meleagris_gallopavo","Turkey"             },
        { "Oryctolagus_cuniculus","Rabbit"            },
        { "Camelus_dromedarius","Camel"              },
        { NULL, NULL }
    };
    for (int i = 0; names[i].latin; i++) {
        if (strcmp(species_id, names[i].latin) == 0)
            return names[i].common;
    }
    return species_id;  /* fallback to Latin name */
}

/* Friendly verdict text */
static const char *friendly_verdict(verdict_t v) {
    switch (v) {
        case PASS: return "HALAL - No haram content detected";
        case FAIL: return "NOT HALAL - Haram content detected";
        default:   return "INCONCLUSIVE - Unable to determine";
    }
}

/* Friendly status text */
static const char *friendly_status(halal_status_t s) {
    switch (s) {
        case HALAL:    return "Halal";
        case HARAM:    return "Haram";
        case MASHBOOH: return "Doubtful";
        default:       return "Unknown";
    }
}

/* Confidence level from cross-marker agreement */
static const char *confidence_label(double agreement) {
    if (agreement >= 0.95) return "Very High";
    if (agreement >= 0.85) return "High";
    if (agreement >= 0.70) return "Moderate";
    return "Low";
}

/* ================================================================== */
/* File dialog                                                         */
/* ================================================================== */
static void open_file_dialog(char *out, size_t out_sz) {
    const char *filters[] = { "*.fq", "*.fastq", "*.fq.gz", "*.fastq.gz",
                               "*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz" };
    const char *result = tinyfd_openFileDialog(
        "Select DNA sample file",              /* title       */
        "",                                    /* default path*/
        8, filters,                            /* filter      */
        "DNA sample files (*.fq *.fa *.gz)",  /* description */
        0                                      /* multi-select*/
    );
    if (result) {
        snprintf(out, out_sz, "%s", result);
    }
}

/* ================================================================== */
/* Locate default index file                                           */
/* ================================================================== */
static void find_default_index(char *path, size_t sz) {
#ifdef __APPLE__
    {
        char buf[1024];
        const char *base = SDL_GetBasePath();
        if (base) {
            snprintf(buf, sizeof(buf), "%s../Resources/default.idx", base);
            FILE *f = fopen(buf, "rb");
            if (f) { fclose(f); snprintf(path, sz, "%s", buf); return; }
        }
    }
#endif
    {
        const char *base = SDL_GetBasePath();
        if (base) {
            char buf[1024];
            snprintf(buf, sizeof(buf), "%sdefault.idx", base);
            FILE *f = fopen(buf, "rb");
            if (f) { fclose(f); snprintf(path, sz, "%s", buf); return; }
        }
    }
    {
        FILE *f = fopen("halal.idx", "rb");
        if (f) { fclose(f); snprintf(path, sz, "halal.idx"); return; }
    }
    path[0] = '\0';
}

/* ================================================================== */
/* Format file size                                                    */
/* ================================================================== */
static void format_file_size(const char *path, char *out, size_t sz) {
    FILE *f = fopen(path, "rb");
    if (!f) { snprintf(out, sz, "-"); return; }
    fseek(f, 0, SEEK_END);
    long bytes = ftell(f);
    fclose(f);
    if (bytes < 1024)
        snprintf(out, sz, "%ld B", bytes);
    else if (bytes < 1024 * 1024)
        snprintf(out, sz, "%.1f KB", bytes / 1024.0);
    else
        snprintf(out, sz, "%.1f MB", bytes / (1024.0 * 1024.0));
}

/* ================================================================== */
/* Draw stacked species bar (bottom of left panel)                     */
/* ================================================================== */
static void draw_species_bar(struct nk_context *ctx, const halal_report_t *r) {
    struct nk_command_buffer *canvas = nk_window_get_canvas(ctx);
    struct nk_rect bounds;
    nk_layout_row_dynamic(ctx, ROW_BAR, 1);
    nk_widget(&bounds, ctx);

    float x = bounds.x;
    for (int i = 0; i < r->n_species; i++) {
        float w = (float)(r->species[i].weight_pct / 100.0) * bounds.w;
        if (w < 1.0f) continue;
        struct nk_color c = status_color(r->species[i].halal_status);
        nk_fill_rect(canvas, nk_rect(x, bounds.y, w, bounds.h), 0, c);
        x += w;
    }
}

/* ================================================================== */
/* Draw the full GUI layout                                            */
/* ================================================================== */
static void draw_gui(struct nk_context *ctx,
                     analysis_context_t *analysis,
                     char *selected_file,
                     char *index_path,
                     int win_w, int win_h)
{
    if (!nk_begin(ctx, "HalalSeq",
                  nk_rect(0, 0, (float)win_w, (float)win_h),
                  NK_WINDOW_NO_SCROLLBAR | NK_WINDOW_BACKGROUND)) {
        nk_end(ctx);
        return;
    }

    /* Title bar */
    nk_layout_row_dynamic(ctx, ROW_TITLE, 1);
    nk_label(ctx, "HalalSeq - Halal Food DNA Authentication", NK_TEXT_CENTERED);

    nk_layout_row_dynamic(ctx, 2, 1);
    nk_spacing(ctx, 1);

    /* Two-column layout */
    float col_widths[] = { (float)LEFT_W, (float)(win_w - LEFT_W - 20) };
    nk_layout_row(ctx, NK_STATIC, (float)(win_h - ROW_TITLE - 24), 2, col_widths);

    /* ============================================================== */
    /* LEFT PANEL                                                      */
    /* ============================================================== */
    if (nk_group_begin(ctx, "input_panel", NK_WINDOW_BORDER)) {

        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        nk_label(ctx, "DNA Sample File:", NK_TEXT_LEFT);

        /* Selected file path */
        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        if (selected_file[0]) {
            const char *fname = strrchr(selected_file, '/');
            if (!fname) fname = strrchr(selected_file, '\\');
            fname = fname ? fname + 1 : selected_file;
            nk_label(ctx, fname, NK_TEXT_LEFT);
        } else {
            nk_label(ctx, "(drop file or click Choose File)",
                     NK_TEXT_LEFT);
        }

        /* File size */
        if (selected_file[0]) {
            char sz_buf[64];
            format_file_size(selected_file, sz_buf, sizeof(sz_buf));
            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            char info[128];
            snprintf(info, sizeof(info), "Size: %s", sz_buf);
            nk_label(ctx, info, NK_TEXT_LEFT);
        }

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        /* Browse button */
        nk_layout_row_dynamic(ctx, ROW_BTN, 1);
        if (nk_button_label(ctx, "Choose File...")) {
            open_file_dialog(selected_file, 1024);
        }

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        /* Analyze button */
        nk_layout_row_dynamic(ctx, ROW_BTN, 1);
        {
            int can_run = (selected_file[0] &&
                          (analysis->state == ANALYSIS_IDLE ||
                           analysis->state == ANALYSIS_DONE ||
                           analysis->state == ANALYSIS_ERROR) &&
                          index_path[0]);
            if (can_run) {
                struct nk_style_button green = ctx->style.button;
                green.normal  = nk_style_item_color(nk_rgb(34, 139, 34));
                green.hover   = nk_style_item_color(nk_rgb(0, 180, 0));
                green.active  = nk_style_item_color(nk_rgb(0, 140, 0));
                green.text_normal  = nk_rgb(255, 255, 255);
                green.text_hover   = nk_rgb(255, 255, 255);
                green.text_active  = nk_rgb(255, 255, 255);
                if (nk_button_label_styled(ctx, &green, "Run Analysis")) {
                    snprintf(analysis->fastq_path,
                             sizeof(analysis->fastq_path), "%s",
                             selected_file);
                    snprintf(analysis->index_path,
                             sizeof(analysis->index_path), "%s",
                             index_path);
                    analysis_start(analysis);
                }
            } else {
                struct nk_style_button grey = ctx->style.button;
                grey.normal = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.hover  = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.active = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.text_normal = nk_rgb(140, 140, 140);
                grey.text_hover  = nk_rgb(140, 140, 140);
                grey.text_active = nk_rgb(140, 140, 140);
                nk_button_label_styled(ctx, &grey, "Run Analysis");
            }
        }

        /* Status line — friendly labels */
        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        {
            const char *status_text;
            switch (analysis->state) {
                case ANALYSIS_IDLE:              status_text = "Ready"; break;
                case ANALYSIS_LOADING_INDEX:     status_text = "Preparing..."; break;
                case ANALYSIS_READING_FASTQ:     status_text = "Reading sample..."; break;
                case ANALYSIS_CLASSIFYING:       status_text = "Identifying species..."; break;
                case ANALYSIS_RUNNING_EM:        status_text = "Calculating amounts..."; break;
                case ANALYSIS_GENERATING_REPORT: status_text = "Generating report..."; break;
                case ANALYSIS_DONE:              status_text = "Analysis complete"; break;
                case ANALYSIS_ERROR:             status_text = "Error occurred"; break;
                default:                         status_text = ""; break;
            }
            nk_label(ctx, status_text, NK_TEXT_LEFT);
        }

        /* Progress bar */
        if (analysis->state > ANALYSIS_IDLE &&
            analysis->state < ANALYSIS_DONE) {
            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            nk_size pv = 0;
            switch (analysis->state) {
                case ANALYSIS_LOADING_INDEX:     pv = 10; break;
                case ANALYSIS_READING_FASTQ:     pv = 25; break;
                case ANALYSIS_CLASSIFYING:       pv = 50; break;
                case ANALYSIS_RUNNING_EM:        pv = 75; break;
                case ANALYSIS_GENERATING_REPORT: pv = 90; break;
                default: break;
            }
            nk_progress(ctx, &pv, 100, NK_FIXED);
        }

        /* Stacked species bar + legend (after results) */
        if (analysis->report && analysis->state == ANALYSIS_DONE) {
            nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
            nk_spacing(ctx, 1);
            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            nk_label(ctx, "Sample Composition:", NK_TEXT_LEFT);
            draw_species_bar(ctx, analysis->report);

            /* Legend: show species names next to color blocks */
            nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
            nk_spacing(ctx, 1);
            for (int i = 0; i < analysis->report->n_species; i++) {
                species_report_t *sp = &analysis->report->species[i];
                if (sp->weight_pct < 0.5) continue;
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char legend[128];
                snprintf(legend, sizeof(legend), "  %s: %.1f%%",
                         friendly_species_name(sp->species_id),
                         sp->weight_pct);
                nk_label_colored(ctx, legend, NK_TEXT_LEFT,
                                 status_color(sp->halal_status));
            }
        }

        nk_group_end(ctx);
    }

    /* ============================================================== */
    /* RIGHT PANEL                                                     */
    /* ============================================================== */
    if (nk_group_begin(ctx, "results_panel", NK_WINDOW_BORDER)) {

        if (analysis->report && analysis->state == ANALYSIS_DONE) {
            halal_report_t *rpt = analysis->report;

            /* Verdict — large, descriptive */
            nk_layout_row_dynamic(ctx, 44, 1);
            {
                struct nk_color vc = verdict_color(rpt->verdict);
                nk_label_colored(ctx, friendly_verdict(rpt->verdict),
                                 NK_TEXT_CENTERED, vc);
            }

            nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
            nk_spacing(ctx, 1);

            /* Column headers — friendly names, 4 columns (dropped Reads%) */
            nk_layout_row_dynamic(ctx, ROW_TABLE, 4);
            nk_label(ctx, "Animal",  NK_TEXT_LEFT);
            nk_label(ctx, "Status",  NK_TEXT_CENTERED);
            nk_label(ctx, "Amount",  NK_TEXT_RIGHT);
            nk_label(ctx, "Range",   NK_TEXT_RIGHT);

            /* Species rows */
            for (int s = 0; s < rpt->n_species; s++) {
                species_report_t *sp = &rpt->species[s];
                if (sp->weight_pct < 0.01 && sp->read_pct < 0.01)
                    continue;

                nk_layout_row_dynamic(ctx, ROW_TABLE, 4);

                /* Common name */
                nk_label(ctx, friendly_species_name(sp->species_id),
                         NK_TEXT_LEFT);

                /* Status (colour-coded) */
                nk_label_colored(ctx, friendly_status(sp->halal_status),
                                 NK_TEXT_CENTERED,
                                 status_color(sp->halal_status));

                /* Amount */
                char buf[64];
                snprintf(buf, sizeof(buf), "%.1f%%", sp->weight_pct);
                nk_label(ctx, buf, NK_TEXT_RIGHT);

                /* Range */
                snprintf(buf, sizeof(buf), "%.1f-%.1f%%",
                         sp->ci_lo, sp->ci_hi);
                nk_label(ctx, buf, NK_TEXT_RIGHT);
            }

            /* Summary — plain English */
            nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
            nk_spacing(ctx, 1);

            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            {
                char meta[256];
                snprintf(meta, sizeof(meta),
                         "%d DNA fragments analyzed",
                         rpt->total_reads);
                nk_label(ctx, meta, NK_TEXT_LEFT);
            }

            if (rpt->cross_marker_agreement > 0) {
                nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                char conf[128];
                snprintf(conf, sizeof(conf), "Confidence: %s",
                         confidence_label(rpt->cross_marker_agreement));
                nk_label(ctx, conf, NK_TEXT_LEFT);
            }

        } else if (analysis->state == ANALYSIS_IDLE) {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_wrap(ctx,
                "Choose a DNA sample file and click "
                "Run Analysis to check halal status.");
        } else if (analysis->state == ANALYSIS_ERROR) {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_colored_wrap(ctx, analysis->error_msg,
                                  nk_rgb(220, 60, 60));
        } else {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_wrap(ctx, "Analyzing your sample, please wait...");
        }

        nk_group_end(ctx);
    }

    nk_end(ctx);
}

/* ================================================================== */
/* Main                                                                */
/* ================================================================== */
int main(int argc, char *argv[]) {
    (void)argc; (void)argv;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_SetHint(SDL_HINT_VIDEO_HIGHDPI_DISABLED, "0");

    SDL_Window *window = SDL_CreateWindow(
        "HalalSeq - Halal Food DNA Authentication",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WINDOW_W, WINDOW_H,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE |
        SDL_WINDOW_ALLOW_HIGHDPI);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(
        window, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    /* --- Compute DPI scale and set renderer scale -------------------- */
    float dpi_scale = 1.0f;
    {
        int render_w, window_w;
        SDL_GetRendererOutputSize(renderer, &render_w, NULL);
        SDL_GetWindowSize(window, &window_w, NULL);
        if (window_w > 0)
            dpi_scale = (float)render_w / (float)window_w;
        if (dpi_scale < 1.0f) dpi_scale = 1.0f;
    }
    SDL_RenderSetScale(renderer, dpi_scale, dpi_scale);

    /* --- Nuklear init ---------------------------------------------- */
    struct nk_context *ctx = nk_sdl_init(window, renderer);
    {
        struct nk_font_atlas *atlas;
        struct nk_font_config cfg = nk_font_config(0);
        struct nk_font *font = NULL;

        cfg.oversample_h = 3;
        cfg.oversample_v = 2;

        float font_size = 18.0f * dpi_scale;

        nk_sdl_font_stash_begin(&atlas);

        static const char *font_paths[] = {
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/Helvetica.ttc",
            "C:\\Windows\\Fonts\\arial.ttf",
            "C:\\Windows\\Fonts\\segoeui.ttf",
            NULL
        };
        for (const char **p = font_paths; *p; p++) {
            font = nk_font_atlas_add_from_file(atlas, *p, font_size, &cfg);
            if (font) break;
        }
        if (!font)
            font = nk_font_atlas_add_default(atlas, font_size, &cfg);

        nk_sdl_font_stash_end();

        font->handle.height /= dpi_scale;
        nk_style_set_font(ctx, &font->handle);
    }

    /* Dark theme */
    {
        struct nk_color table[NK_COLOR_COUNT];
        table[NK_COLOR_TEXT]                    = nk_rgba(210, 210, 210, 255);
        table[NK_COLOR_WINDOW]                  = nk_rgba(35, 35, 38, 255);
        table[NK_COLOR_HEADER]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_BORDER]                  = nk_rgba(65, 65, 70, 255);
        table[NK_COLOR_BUTTON]                  = nk_rgba(60, 60, 65, 255);
        table[NK_COLOR_BUTTON_HOVER]            = nk_rgba(75, 75, 80, 255);
        table[NK_COLOR_BUTTON_ACTIVE]           = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_TOGGLE]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_TOGGLE_HOVER]            = nk_rgba(55, 55, 60, 255);
        table[NK_COLOR_TOGGLE_CURSOR]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SELECT]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_SELECT_ACTIVE]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SLIDER]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_SLIDER_CURSOR]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SLIDER_CURSOR_HOVER]     = nk_rgba(60, 180, 60, 255);
        table[NK_COLOR_SLIDER_CURSOR_ACTIVE]    = nk_rgba(34, 140, 34, 255);
        table[NK_COLOR_PROPERTY]                = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_EDIT]                    = nk_rgba(45, 45, 50, 255);
        table[NK_COLOR_EDIT_CURSOR]             = nk_rgba(210, 210, 210, 255);
        table[NK_COLOR_COMBO]                   = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_CHART]                   = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_CHART_COLOR]             = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_CHART_COLOR_HIGHLIGHT]   = nk_rgba(255, 0, 0, 255);
        table[NK_COLOR_SCROLLBAR]               = nk_rgba(40, 40, 45, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR]        = nk_rgba(60, 60, 65, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR_HOVER]  = nk_rgba(75, 75, 80, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR_ACTIVE] = nk_rgba(55, 55, 60, 255);
        table[NK_COLOR_TAB_HEADER]              = nk_rgba(50, 50, 55, 255);
        nk_style_from_table(ctx, table);
    }

    /* --- App state ------------------------------------------------- */
    char selected_file[1024] = {0};
    char index_path[1024] = {0};

    analysis_context_t analysis;
    analysis_init(&analysis);

    find_default_index(index_path, sizeof(index_path));

    /* --- Main loop ------------------------------------------------- */
    int running = 1;
    while (running) {
        SDL_Event evt;
        nk_input_begin(ctx);
        while (SDL_PollEvent(&evt)) {
            if (evt.type == SDL_QUIT) {
                running = 0;
            }
            if (evt.type == SDL_DROPFILE) {
                snprintf(selected_file, sizeof(selected_file),
                         "%s", evt.drop.file);
                SDL_free(evt.drop.file);
            }
            nk_sdl_handle_event(&evt);
        }
        nk_input_end(ctx);

        int win_w, win_h;
        SDL_GetWindowSize(window, &win_w, &win_h);

        draw_gui(ctx, &analysis, selected_file, index_path, win_w, win_h);

        SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);
        SDL_RenderClear(renderer);
        nk_sdl_render(NK_ANTI_ALIASING_ON);
        SDL_RenderPresent(renderer);
    }

    analysis_cleanup(&analysis);
    nk_sdl_shutdown();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
