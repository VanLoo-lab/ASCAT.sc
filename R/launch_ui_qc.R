launch_ui_qc <- function(res)
{
    if (!requireNamespace("bslib", quietly = TRUE))
        stop("Package 'bslib' is required: install.packages('bslib')")

    library(shiny)
    library(bslib)
    library(plotly)
    library(DT)

    # ── Data preparation ─────────────────────────────────────────
    if (is.null(res) || is.null(res$QC_metrics))
        stop("res$QC_metrics not found. Run res <- getQCs(res) first.")

    qc_df <- res$QC_metrics
    rn <- rownames(qc_df)
    qc_df$cell_id <- if (!is.null(rn) &&
                         !identical(rn, as.character(seq_len(nrow(qc_df))))) rn
                      else paste0("cell_", seq_len(nrow(qc_df)))

    qc_df <- qc_df[, colSums(!is.na(qc_df)) > 0, drop = FALSE]
    if (!"cell_id" %in% colnames(qc_df))
        qc_df$cell_id <- rownames(qc_df)
    rownames(qc_df) <- NULL

    numeric_cols <- names(qc_df)[vapply(qc_df, is.numeric, logical(1))]

    # ── Human-readable names & descriptions ──────────────────────
    metric_labels <- c(
        distance_integer_logR = "LogR Integer Distance",
        distance_integer_BAF  = "BAF Integer Distance",
        MAPD_gc_corrected     = "MAPD (GC-corrected)",
        spikiness             = "Spikiness",
        MBRSM_dispersion      = "Segment Dispersion",
        breakpoints           = "Breakpoints",
        state_mode            = "Modal CN State",
        autocorr              = "Autocorrelation",
        MSRSI_non_integerness = "Non-integerness (MSRSI)"
    )

    metric_help <- c(
        distance_integer_logR = "Weighted distance of logR copy numbers from nearest integer. Lower = better fit.",
        distance_integer_BAF  = "Weighted distance of BAF allelic states from integers. Lower = cleaner calls.",
        MAPD_gc_corrected     = "Median Absolute Pairwise Difference of GC-corrected logR. Measures bin-level noise.",
        spikiness             = "Second-order difference of logR. Detects isolated outlier bins. Lower = smoother.",
        MBRSM_dispersion      = "Scaled MAD of within-segment residuals. Measures intra-segment noise.",
        breakpoints           = "Number of copy-number state transitions across the genome.",
        state_mode            = "Most frequent copy number state. Typically 2 for near-diploid cells.",
        autocorr              = "Lag-1 autocorrelation of segment residuals. High values suggest systematic bias.",
        MSRSI_non_integerness = "Scaled MAD of segment deviations from integer-expected logR. Lower = more integer-like."
    )

    display_name <- function(col) {
        lbl <- metric_labels[col]
        ifelse(is.na(lbl), gsub("_", " ", col), lbl)
    }

    axis_choices <- stats::setNames(numeric_cols, display_name(numeric_cols))

    # ── Slider precision (2-3 significant digits) ────────────────
    slider_params <- function(rng) {
        span <- diff(rng)
        if (is.na(span) || span == 0) {
            mid <- if (is.na(rng[1])) 0 else rng[1]
            return(list(min = mid - 1, max = mid + 1, step = 0.01))
        }
        step <- signif(span / 100, 1)
        lo   <- floor(rng[1] / step) * step
        hi   <- ceiling(rng[2] / step) * step
        list(min = lo, max = hi, step = step)
    }

    # ── Default thresholds (median ± 3*MAD) ────────────────────
    default_bounds <- lapply(
        stats::setNames(numeric_cols, numeric_cols),
        function(col) {
            x   <- qc_df[[col]]
            x   <- x[!is.na(x)]
            rng <- range(x)
            p   <- slider_params(rng)
            if (length(x) < 3) return(c(p$min, p$max))
            med <- median(x)
            mad_val <- stats::mad(x)
            if (mad_val == 0) return(c(p$min, p$max))
            lo <- max(p$min, floor((med - 3 * mad_val) / p$step) * p$step)
            hi <- min(p$max, ceiling((med + 3 * mad_val) / p$step) * p$step)
            c(lo, hi)
        }
    )

    # ── Colours (tuned for dark backgrounds) ─────────────────────
    col_pass       <- "#5dade2"
    col_fail       <- "#ec7063"
    col_pass_alpha <- "rgba(93,173,226,0.55)"
    col_fail_alpha <- "rgba(236,112,99,0.55)"
    col_median     <- "#f39c12"
    col_text       <- "#dee2e6"
    col_grid       <- "rgba(255,255,255,0.08)"

    # ── Custom CSS (dark mode) ───────────────────────────────────
    app_css <- "
        /* ── Dark scrollbars (native + fallback) ── */
        html, body { color-scheme: dark; }
        ::-webkit-scrollbar { width: 7px; height: 7px; }
        ::-webkit-scrollbar-track { background: #2b2b2b; }
        ::-webkit-scrollbar-thumb { background: #555; border-radius: 4px; }
        ::-webkit-scrollbar-thumb:hover { background: #777; }

        .value-box .value-box-value { font-variant-numeric: tabular-nums; }
        .sidebar .form-group { margin-bottom: 0.6rem; }
        .sidebar .control-label { font-size: 0.82rem; }
        .section-label {
            text-transform: uppercase; font-weight: 700; font-size: 0.7rem;
            letter-spacing: 1.2px; color: #95a5a6; margin-bottom: 0.75rem;
        }
        .metric-help {
            cursor: help; opacity: 0.45; margin-left: 4px; font-size: 0.75rem;
        }
        .metric-help:hover { opacity: 0.9; }
        .card, .bslib-grid { margin-bottom: 1.3rem; }
        .card { border-color: rgba(255,255,255,0.08); }
        .card-header { border-bottom-color: rgba(255,255,255,0.08); }
        hr { border-color: rgba(255,255,255,0.1); }

        /* ── Slider theming (ion.rangeSlider) ── */
        .irs--shiny .irs-bar {
            background: #5dade2; border-top: 1px solid #5dade2;
            border-bottom: 1px solid #5dade2;
        }
        .irs--shiny .irs-handle {
            background: #dee2e6; border: 2px solid #5dade2;
            box-shadow: none;
        }
        .irs--shiny .irs-handle:hover { background: #fff; }
        .irs--shiny .irs-line {
            background: #3a3f47; border-color: #3a3f47;
        }
        .irs--shiny .irs-from,
        .irs--shiny .irs-to,
        .irs--shiny .irs-single {
            background: #5dade2; color: #fff; font-size: 0.75rem;
            padding: 1px 5px;
        }
        .irs--shiny .irs-min,
        .irs--shiny .irs-max {
            background: #3a3f47; color: #7f8c8d; font-size: 0.7rem;
        }
        .irs--shiny .irs-grid-text { color: #6c757d; font-size: 0.65rem; }
        .irs--shiny .irs-grid-pol { background: #4a4a4a; }
    "

    # ── UI ───────────────────────────────────────────────────────
    ui <- page_sidebar(
        title = div(
            style = "display:flex; align-items:center; gap:10px;",
            span(style = "font-weight:700; letter-spacing:0.5px;", "ASCAT.sc"),
            span(style = "opacity:0.35; font-weight:300;", "|"),
            span(style = "font-weight:400; opacity:0.6;",
                 "Quality Control Dashboard")
        ),
        fillable = FALSE,
        theme = bs_theme(
            version    = 5,
            bootswatch = "darkly",
            primary    = "#5dade2",
            success    = "#2ecc71",
            danger     = "#e74c3c",
            info       = "#3498db"
        ),
        tags$head(tags$style(HTML(app_css))),

        sidebar = sidebar(
            width = 340,

            # ---- Axes ----
            tags$div(class = "section-label", "Plot Axes"),
            selectInput("xvar", "X axis", axis_choices,
                        selected = numeric_cols[min(1, length(numeric_cols))]),
            selectInput("yvar", "Y axis", axis_choices,
                        selected = numeric_cols[min(2, length(numeric_cols))]),
            selectInput("zvar", "Z axis", axis_choices,
                        selected = numeric_cols[min(3, length(numeric_cols))]),

            hr(style = "margin:0.8rem 0;"),

            # ---- Gates ----
            tags$div(class = "section-label", "Gate Thresholds"),
            tags$div(
                class = "gate-scroll",
                style = "max-height:380px; overflow-y:auto; padding:0 8px 0 14px;",
                uiOutput("gate_sliders")
            ),
            actionButton("reset_gates", "Reset All Gates",
                         class = "btn-outline-light btn-sm w-100 mt-2"),

            hr(style = "margin:0.8rem 0;"),

            # ---- Export ----
            tags$div(class = "section-label", "Export"),
            downloadButton("download_pass",  "Passing Cells",
                           class = "btn-sm btn-outline-success w-100 mb-2"),
            downloadButton("download_fail",  "Failing Cells",
                           class = "btn-sm btn-outline-danger  w-100 mb-2"),
            downloadButton("download_bounds", "Gate Bounds",
                           class = "btn-sm btn-outline-light w-100")
        ),

        # ── Main area ────────────────────────────────────────────
        layout_columns(
            col_widths = c(4, 4, 4),
            value_box("Total Cells", textOutput("n_total"),
                      theme = "secondary"),
            value_box("Passing",     textOutput("n_pass"),
                      theme = "success"),
            value_box("Failing",     textOutput("n_fail"),
                      theme = "danger")
        ),

        card(
            card_header(
                class = "d-flex justify-content-between align-items-center",
                tags$span(class = "fw-bold", "3D QC Explorer"),
                tags$small(
                    class = "text-muted",
                    "Drag to rotate • Scroll to zoom • Right-drag to pan"
                )
            ),
            plotlyOutput("plot3d", height = "700px")
        ),

        layout_columns(
            col_widths = c(4, 4, 4),
            card(
                card_header(class = "py-2 fw-semibold",
                            textOutput("hist_title_x")),
                card_body(class = "p-1",
                          plotlyOutput("hist_x", height = "250px"))
            ),
            card(
                card_header(class = "py-2 fw-semibold",
                            textOutput("hist_title_y")),
                card_body(class = "p-1",
                          plotlyOutput("hist_y", height = "250px"))
            ),
            card(
                card_header(class = "py-2 fw-semibold",
                            textOutput("hist_title_z")),
                card_body(class = "p-1",
                          plotlyOutput("hist_z", height = "250px"))
            )
        ),

        card(
            card_header(
                class = "d-flex justify-content-between align-items-center",
                tags$span(class = "fw-bold", "Cell Data"),
                tags$small(class = "text-muted",
                           textOutput("table_summary", inline = TRUE))
            ),
            DTOutput("table")
        )
    )

    # ── Server ───────────────────────────────────────────────────
    server <- function(input, output, session)
    {
        bounds <- reactiveValues()

        observe({
            for (col in numeric_cols) {
                if (is.null(bounds[[col]]))
                    bounds[[col]] <- default_bounds[[col]]
            }
        })

        # ---- Gate sliders ----
        output$gate_sliders <- renderUI({
            sliders <- lapply(numeric_cols, function(col) {
                rng  <- range(qc_df[[col]], na.rm = TRUE)
                p    <- slider_params(rng)
                curr <- default_bounds[[col]]

                help <- metric_help[col]
                lbl  <- if (!is.na(help)) {
                    tags$span(
                        display_name(col),
                        tags$span(class = "metric-help", title = help,
                                  "ⓘ")
                    )
                } else {
                    display_name(col)
                }

                sliderInput(
                    paste0("gate_", col),
                    label = lbl,
                    min   = p$min,
                    max   = p$max,
                    value = curr,
                    step  = p$step
                )
            })
            tagList(sliders)
        })

        observe({
            for (col in numeric_cols) {
                val <- input[[paste0("gate_", col)]]
                if (!is.null(val)) bounds[[col]] <- val
            }
        })

        observeEvent(input$reset_gates, {
            for (col in numeric_cols) {
                bounds[[col]] <- default_bounds[[col]]
                updateSliderInput(session, paste0("gate_", col),
                                  value = default_bounds[[col]])
            }
        })

        # ---- Classification ----
        classified <- reactive({
            df <- qc_df
            df$pass <- TRUE
            for (col in numeric_cols) {
                b <- bounds[[col]]
                if (!is.null(b)) {
                    v <- df[[col]]
                    df$pass <- df$pass & (is.na(v) | (v >= b[1] & v <= b[2]))
                }
            }
            df$status <- ifelse(df$pass, "Pass", "Fail")
            df
        })

        # ---- Summary ----
        output$n_total <- renderText(nrow(qc_df))
        output$n_pass  <- renderText({
            n <- sum(classified()$pass)
            paste0(n, " (", round(100 * n / nrow(qc_df), 1), "%)")
        })
        output$n_fail  <- renderText({
            n <- sum(!classified()$pass)
            paste0(n, " (", round(100 * n / nrow(qc_df), 1), "%)")
        })
        output$table_summary <- renderText({
            df <- classified()
            paste0(sum(df$pass), " passing, ",
                   sum(!df$pass), " failing of ",
                   nrow(df), " cells")
        })

        # ---- 3D plot ----
        output$plot3d <- renderPlotly({
            df  <- classified()
            msz <- if (nrow(df) < 100) 5 else if (nrow(df) < 500) 4 else 3

            pass_df <- df[df$pass,  , drop = FALSE]
            fail_df <- df[!df$pass, , drop = FALSE]

            hover_text <- function(d) {
                paste0("<b>", d$cell_id, "</b><br>",
                       display_name(input$xvar), ": ",
                       signif(d[[input$xvar]], 3), "<br>",
                       display_name(input$yvar), ": ",
                       signif(d[[input$yvar]], 3), "<br>",
                       display_name(input$zvar), ": ",
                       signif(d[[input$zvar]], 3))
            }

            dark_axis <- function(ttl) {
                list(title       = ttl,
                     color       = col_text,
                     gridcolor   = col_grid,
                     zerolinecolor = col_grid,
                     backgroundcolor = "rgba(0,0,0,0)")
            }

            p <- plot_ly()

            if (nrow(pass_df) > 0) {
                p <- p %>% add_trace(
                    x = pass_df[[input$xvar]],
                    y = pass_df[[input$yvar]],
                    z = pass_df[[input$zvar]],
                    type = "scatter3d", mode = "markers",
                    marker = list(size = msz, color = col_pass,
                                  opacity = 0.7,
                                  line = list(width = 0)),
                    text = hover_text(pass_df),
                    hoverinfo = "text", name = "Pass"
                )
            }
            if (nrow(fail_df) > 0) {
                p <- p %>% add_trace(
                    x = fail_df[[input$xvar]],
                    y = fail_df[[input$yvar]],
                    z = fail_df[[input$zvar]],
                    type = "scatter3d", mode = "markers",
                    marker = list(size = msz + 1, color = col_fail,
                                  opacity = 0.9,
                                  line = list(width = 0.5,
                                              color = "#c0392b")),
                    text = hover_text(fail_df),
                    hoverinfo = "text", name = "Fail"
                )
            }

            p %>% layout(
                font  = list(color = col_text),
                scene = list(
                    xaxis  = dark_axis(display_name(input$xvar)),
                    yaxis  = dark_axis(display_name(input$yvar)),
                    zaxis  = dark_axis(display_name(input$zvar)),
                    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
                ),
                legend = list(orientation = "h",
                              x = 0.5, xanchor = "center",
                              y = 1.02, yanchor = "bottom",
                              font = list(color = col_text)),
                margin = list(t = 40, b = 0, l = 0, r = 0),
                paper_bgcolor = "rgba(0,0,0,0)",
                plot_bgcolor  = "rgba(0,0,0,0)"
            ) %>% config(displaylogo = FALSE)
        })

        # ---- Histograms ----
        make_histogram <- function(col) {
            req(col)
            df <- classified()

            pass_vals <- df[[col]][df$pass]
            fail_vals <- df[[col]][!df$pass]
            med_pass  <- median(pass_vals, na.rm = TRUE)

            p <- plot_ly() %>%
                add_histogram(
                    x = pass_vals, name = "Pass",
                    marker = list(color = col_pass_alpha,
                                  line  = list(color = col_pass,
                                               width = 1))
                ) %>%
                add_histogram(
                    x = fail_vals, name = "Fail",
                    marker = list(color = col_fail_alpha,
                                  line  = list(color = col_fail,
                                               width = 1))
                )

            shapes      <- list()
            annotations <- list()

            if (!is.na(med_pass)) {
                shapes <- list(list(
                    type = "line",
                    x0 = med_pass, x1 = med_pass,
                    y0 = 0, y1 = 1, yref = "paper",
                    line = list(color = col_median, width = 2,
                                dash = "dash")
                ))
                annotations <- list(list(
                    x = med_pass, y = 1, yref = "paper",
                    text      = paste0("median: ", signif(med_pass, 3)),
                    showarrow = FALSE,
                    yanchor   = "bottom",
                    font      = list(color = col_median, size = 11)
                ))
            }

            p %>% layout(
                barmode       = "overlay",
                font          = list(color = col_text),
                xaxis         = list(title = display_name(col),
                                     color = col_text,
                                     gridcolor = col_grid,
                                     zeroline  = FALSE),
                yaxis         = list(title = "Count",
                                     color = col_text,
                                     gridcolor = col_grid,
                                     zeroline  = FALSE),
                shapes        = shapes,
                annotations   = annotations,
                showlegend    = FALSE,
                margin        = list(t = 20, b = 40, l = 45, r = 10),
                paper_bgcolor = "rgba(0,0,0,0)",
                plot_bgcolor  = "rgba(0,0,0,0)"
            ) %>% config(displayModeBar = FALSE)
        }

        output$hist_x <- renderPlotly(make_histogram(input$xvar))
        output$hist_y <- renderPlotly(make_histogram(input$yvar))
        output$hist_z <- renderPlotly(make_histogram(input$zvar))

        output$hist_title_x <- renderText(display_name(input$xvar))
        output$hist_title_y <- renderText(display_name(input$yvar))
        output$hist_title_z <- renderText(display_name(input$zvar))

        # ---- Data table ----
        output$table <- renderDT({
            df   <- classified()
            disp <- df[, c("cell_id", "status", numeric_cols), drop = FALSE]
            for (col in numeric_cols)
                disp[[col]] <- signif(disp[[col]], 3)

            datatable(
                disp, rownames = FALSE, filter = "top",
                options = list(
                    pageLength = 15, scrollX = TRUE,
                    dom = "lftip",
                    order = list(list(1, "asc")),
                    initComplete = JS(
                        "function(settings, json) {",
                        "  $(this.api().table().container())",
                        "    .css({'color':'#dee2e6','background':'#303030'});",
                        "  $(this.api().table().header())",
                        "    .css({'background':'#3a3a3a','color':'#dee2e6'});",
                        "}"
                    )
                ),
                class = "compact stripe hover"
            ) %>%
                formatStyle("status",
                            backgroundColor = styleEqual(
                                c("Pass", "Fail"),
                                c("rgba(46,204,113,0.18)",
                                  "rgba(231,76,60,0.18)")),
                            color = styleEqual(
                                c("Pass", "Fail"),
                                c("#2ecc71", "#e74c3c")),
                            fontWeight = "bold")
        })

        # ---- Downloads ----
        output$download_pass <- downloadHandler(
            filename = function() paste0("qc_pass_", Sys.Date(), ".csv"),
            content  = function(file) {
                df <- classified()
                out <- df[df$pass, !names(df) %in% c("pass", "status"),
                          drop = FALSE]
                write.csv(out, file, row.names = FALSE)
            }
        )
        output$download_fail <- downloadHandler(
            filename = function() paste0("qc_fail_", Sys.Date(), ".csv"),
            content  = function(file) {
                df <- classified()
                out <- df[!df$pass, !names(df) %in% c("pass", "status"),
                          drop = FALSE]
                write.csv(out, file, row.names = FALSE)
            }
        )
        output$download_bounds <- downloadHandler(
            filename = function() paste0("qc_bounds_", Sys.Date(), ".csv"),
            content  = function(file) {
                bounds_df <- data.frame(
                    metric       = numeric_cols,
                    display_name = display_name(numeric_cols),
                    gate_min = vapply(numeric_cols,
                                     function(x) bounds[[x]][1],
                                     numeric(1)),
                    gate_max = vapply(numeric_cols,
                                     function(x) bounds[[x]][2],
                                     numeric(1)),
                    data_min = vapply(numeric_cols,
                                     function(x) min(qc_df[[x]],
                                                     na.rm = TRUE),
                                     numeric(1)),
                    data_max = vapply(numeric_cols,
                                     function(x) max(qc_df[[x]],
                                                     na.rm = TRUE),
                                     numeric(1)),
                    stringsAsFactors = FALSE
                )
                write.csv(bounds_df, file, row.names = FALSE)
            }
        )
    }

    shinyApp(ui, server)
}
