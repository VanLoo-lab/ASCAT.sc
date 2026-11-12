# =============================================================================
# launch_ui_qc(): interactive quality control of ASCAT.sc single-cell results
# -----------------------------------------------------------------------------
#   qc <- getQCs(res)
#   launch_ui_qc(qc)       # tune the thresholds, then click "Done"; qc is updated:
#   qc$pass                # ids of the cells that pass
#   qc$fail                # ids of the cells that fail
#   qc$cells               # one row per cell: status, failed metrics, values
#   qc$QC_session$thresholds  # thresholds used, in original units
#   launch_ui_qc(qc)       # re-opens with the same thresholds and decisions
#
# The results are added to the object you pass in (nothing is removed), and
# it is also returned, so qc <- launch_ui_qc(qc) works the same way.
# =============================================================================

launch_ui_qc <- function(res, theme = c("light", "dark"), previous = res$QC_session, ...)
{
  # Name of the object passed in, so it can be updated where it lives
  caller   <- parent.frame()
  res_name <- if (is.name(substitute(res))) as.character(substitute(res)) else NULL
  theme <- match.arg(theme)
  needed <- c("bslib", "plotly", "DT", "htmltools")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("launch_ui_qc() needs ", paste(missing, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0("\"", missing, "\"", collapse = ", "), "))",
         call. = FALSE)
  if (is.null(res) || is.null(res$QC_metrics))
    stop("res$QC_metrics not found. Run res <- getQCs(res) first.", call. = FALSE)
  # previous can also be an object from an earlier launch_ui_qc() call
  if (!is.null(previous$QC_session)) previous <- previous$QC_session

  qd <- .qc_prepare(res, previous)

  # Shared across browser sessions: a page refresh keeps the gates, and the
  # final state is what launch_ui_qc() returns.
  st <- new.env(parent = emptyenv())
  st$defaults <- .qc_initial_state(qd)
  st$state    <- .qc_initial_state(qd, previous = previous)

  label <- if (is.null(res_name)) "res" else res_name
  app <- shiny::shinyApp(
    ui     = function(req) .qc_ui(qd, st$state, theme),
    server = .qc_server(qd, res, st, theme, label)
  )
  tryCatch(shiny::runApp(app, ...), interrupt = function(e) invisible(NULL))

  res <- .qc_attach(res, .qc_result(qd, st$state))
  # Update the caller's object as well, so launch_ui_qc(qc) needs no assignment
  saved <- !is.null(res_name) && exists(res_name, envir = caller, inherits = FALSE) &&
    !inherits(try(assign(res_name, res, envir = caller), silent = TRUE), "try-error")
  message(sprintf("ASCAT.sc QC: %s of %s cells pass (%s).", format(length(res$pass), big.mark = ","),
                  format(qd$n, big.mark = ","), .qc_pct(length(res$pass), qd$n)))
  message(if (saved)
    sprintf("Saved in %1$s: %1$s$pass, %1$s$fail and %1$s$cells. Run launch_ui_qc(%1$s) again to continue where you left off.",
            res_name)
  else "Keep the result with res <- launch_ui_qc(res), then use res$pass, res$fail and res$cells.")
  invisible(res)
}

# ── UI ───────────────────────────────────────────────────────────────────────
# A segmented control that is a regular Shiny radio input.
.qc_segmented <- function(id, choices, selected = choices[[1]], label = NULL)
{
  shiny::div(
    id = id, class = "shiny-input-radiogroup shiny-input-container qc-seg",
    role = "radiogroup", `aria-label` = label,
    lapply(seq_along(choices), function(i) {
      v <- choices[[i]]
      shiny::tags$label(
        shiny::tags$input(type = "radio", name = id, value = v,
                          checked = if (identical(v, selected)) NA else NULL),
        shiny::span(names(choices)[i]))
    }))
}

.qc_ui <- function(qd, init, theme)
{
  tags <- shiny::tags; div <- shiny::div; span <- shiny::span; icon <- shiny::icon
  M   <- qd$metrics
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

  grouped <- lapply(split(M, factor(M$group, levels = unique(M$group))),
                    function(g) stats::setNames(g$metric, g$label))
  pick <- function(pref, avoid = character()) {
    hit <- setdiff(intersect(pref, M$metric), avoid)
    if (length(hit)) hit[1] else setdiff(M$metric, avoid)[1]
  }
  x0 <- pick(c("n_reads_total", "MAPD_gc_corrected"))
  y0 <- pick(c("MAPD_gc_corrected", "spikiness"), x0)
  z0 <- pick(c("MSRSI_non_integerness", "breakpoints", "spikiness"), c(x0, y0))
  if (is.na(z0)) z0 <- x0
  has_baf <- "distance_integer_BAF" %in% M$metric

  field <- function(label, ...) div(class = "qc-field", span(class = "qc-label", label), ...)
  card_head <- function(title, sub = NULL, right = NULL)
    div(class = "qc-card-head",
        div(class = "qc-card-head-left", span(class = "qc-card-title", title),
            if (!is.null(sub)) span(class = "qc-card-sub", sub)),
        right)
  legend <- function(...) div(class = "qc-legend", ...)
  tab_title <- function(ic, txt) shiny::tagList(icon(ic), span(class = "qc-tab-text", txt))

  # ---- Header ----
  dl_item <- function(id, ic, txt)
    tags$li(shiny::downloadLink(id, shiny::tagList(icon(ic), txt), class = "dropdown-item"))
  header <- tags$header(
    class = "qc-header",
    div(class = "qc-brand", .qc_logo(),
        div(class = "qc-title", "ASCAT.sc ", span("Quality control"))),
    div(class = "qc-pills",
        span(class = "qc-pill", fmt(qd$n), " cells"),
        span(class = "qc-pill", nrow(M), " metrics"),
        if (is.finite(qd$n_bins)) span(class = "qc-pill", fmt(qd$n_bins), " bins"),
        if (has_baf) span(class = "qc-pill", "Allele-specific")),
    div(class = "qc-header-actions",
        tags$button(type = "button", class = "qc-btn qc-btn-icon", `data-qc-toggle-sidebar` = NA,
                    title = "Show or hide the threshold panel", `aria-pressed` = "false",
                    icon("table-columns")),
        div(class = "dropdown",
            tags$button(type = "button", class = "qc-btn dropdown-toggle",
                        `data-bs-toggle` = "dropdown", `aria-expanded` = "false",
                        icon("arrow-up-from-bracket"), "Export"),
            tags$ul(class = "dropdown-menu dropdown-menu-end",
                    tags$li(tags$h6(class = "dropdown-header", "CSV")),
                    dl_item("dl_pass", "circle-check", "Passing cells"),
                    dl_item("dl_fail", "circle-xmark", "Failing cells"),
                    dl_item("dl_table", "table", "All cells with status"),
                    dl_item("dl_gates", "sliders", "Thresholds"))),
        bslib::input_dark_mode(id = "dark_mode", mode = theme),
        shiny::actionButton("done", "Done", icon = icon("check"), class = "qc-btn qc-btn-primary",
                            title = "Close the app and return the result to R")))

  # ---- Sidebar: gates ----
  gate_ui <- function(m) {
    sp <- M[m, ]
    div(class = "qc-gate", `data-metric` = m,
        div(class = "qc-gate-head",
            bslib::input_switch(paste0("on_", m), NULL, value = isTRUE(init$on[[m]])),
            span(class = "qc-gate-name", sp$label),
            if (!is.na(sp$help))
              bslib::tooltip(span(class = "qc-help", tabindex = "0", icon("circle-info")),
                             paste0(sp$help,
                                    if (sp$n_na) sprintf(" %d cells have no value and always pass.", sp$n_na)),
                             placement = "right")),
        div(class = "qc-gate-body",
            shiny::uiOutput(paste0("gs_", m)),
            shiny::sliderInput(paste0("gate_", m), NULL, min = sp$min, max = sp$max,
                               value = init$bounds[[m]], step = sp$step, ticks = FALSE,
                               width = "100%")))
  }
  groups <- unique(M$group)
  gate_sections <- lapply(groups, function(g) shiny::tagList(
    div(class = "qc-label qc-group-label", g),
    lapply(M$metric[M$group == g], gate_ui)))

  flag_ui <- if (nrow(qd$flags)) div(
    class = "qc-side-section",
    div(class = "qc-label qc-group-label", paste0("Flags ", .qc_s[["mid"]], " exclude cells that are")),
    lapply(qd$flags$flag, function(f) {
      fr <- qd$flags[f, ]
      div(class = "qc-flag",
          bslib::input_switch(paste0("flag_", f), NULL, value = isTRUE(init$flags_on[[f]])),
          span(fr$label),
          if (!is.na(fr$help))
            bslib::tooltip(span(class = "qc-help", tabindex = "0", icon("circle-info")), fr$help,
                           placement = "right"),
          span(class = "qc-flag-n", fmt(fr$n), " cells"))
    }))

  sidebar <- tags$aside(
    class = "qc-sidebar", `aria-label` = "Thresholds",
    div(class = "qc-side-section",
        div(class = "qc-label", "Thresholds"),
        div(class = "qc-auto",
            span(paste0("median ", .qc_s[["pm"]])),
            shiny::numericInput("k", NULL, value = init$k, min = 0.5, max = 10, step = 0.5,
                                width = "56px"),
            span(paste0(.qc_s[["times"]], " MAD")),
            div(class = "qc-auto-actions",
                shiny::actionButton("apply_k", "Apply", class = "qc-btn qc-btn-sm",
                                    title = "Recompute every threshold in use from the median and MAD")))),
    gate_sections,
    flag_ui,
    div(class = "qc-side-foot",
        shiny::actionButton("reset", "Reset thresholds", icon = icon("rotate-left"),
                            class = "qc-btn qc-btn-sm")),
    div(class = "qc-kbd-hints",
        div(tags$kbd(.qc_s[["larr"]]), " ", tags$kbd(.qc_s[["rarr"]]), " previous / next cell"),
        div(tags$kbd("lasso"), " select cells in Explore")))

  # ---- Tabs ----
  overview <- shiny::tagList(
    shiny::uiOutput("kpis"),
    div(class = "qc-grid-2",
        div(class = "qc-card",
            card_head("Pass and fail per metric", "cells within / outside each threshold"),
            shiny::uiOutput("impact")),
        div(class = "qc-card",
            card_head("Metric correlations", paste0("Spearman ", .qc_s[["rho"]], ", strength ", .qc_s[["mid"]], " click a square to plot the pair")),
            plotly::plotlyOutput("corr", height = "340px"))),
    div(class = "qc-card", style = "margin-top:16px",
        card_head("Distributions", "all cells, stacked by overall status",
                  legend(span(tags$i(class = "pass"), "Pass"), span(tags$i(class = "fail"), "Fail"),
                         span(tags$i(class = "line"), "threshold"))),
        div(class = "qc-dists", lapply(M$metric, function(m)
          div(class = "qc-dist",
              div(class = "qc-dist-head", span(class = "qc-dist-name", M[m, "label"]),
                  shiny::textOutput(paste0("distgate_", m), inline = TRUE,
                                    container = function(...) span(class = "qc-dist-gate", ...))),
              plotly::plotlyOutput(paste0("hist_", m), height = "128px"))))))

  axis_select <- function(id, sel)
    div(class = "qc-select-sm", shiny::selectInput(id, NULL, grouped, selected = sel, width = "190px"))
  explore <- shiny::tagList(
    div(class = "qc-toolbar",
        field("X", axis_select("xvar", x0)),
        field("Y", axis_select("yvar", y0)),
        shiny::conditionalPanel("input.dim === '3d'", field("Z", axis_select("zvar", z0))),
        field("Colour", div(class = "qc-select-sm",
                            shiny::selectInput("color_by", NULL,
                                               c(list(Status = c("Pass / fail" = "status")), grouped),
                                               selected = "status", width = "170px"))),
        div(class = "qc-spacer"),
        .qc_segmented("dim", c("3D" = "3d", "2D" = "2d"), "3d", "Plot dimensions")),
    shiny::uiOutput("selbar"),
    div(class = "qc-card",
        plotly::plotlyOutput("explore", height = "600px"),
        shiny::conditionalPanel(
          "input.dim === '3d'",
          div(class = "qc-hint", icon("cube"),
              paste0("Drag to rotate ", .qc_s[["mid"]], " scroll to zoom ", .qc_s[["mid"]], " right-drag to pan ", .qc_s[["mid"]], " click a cell to see its profile."),
              " The box marks the thresholds on the three axes; switch to 2D to select cells with the lasso.")),
        shiny::conditionalPanel(
          "input.dim === '2d'",
          div(class = "qc-hint", icon("hand-pointer"),
              "Click a cell to see its profile. Use the box or lasso tool to select several cells."))),
    shiny::conditionalPanel(
      "output.has_sel",
      div(class = "qc-card",
          div(class = "qc-card-head", shiny::uiOutput("explore_sel_head", class = "qc-explore-sel-head"),
              shiny::actionButton("to_inspector", "Open in Inspector", icon = icon("arrow-right"),
                                  class = "qc-btn qc-btn-sm")),
          div(class = "qc-cn-row", .qc_cn_legend()),
          plotly::plotlyOutput("explore_profile", height = "250px"))))

  cells <- shiny::tagList(
    div(class = "qc-toolbar",
        .qc_segmented("tbl_filter", c("All" = "all", "Passing" = "pass", "Failing" = "fail",
                                      "Manual" = "manual"), "all", "Filter cells"),
        div(class = "qc-spacer"),
        span(class = "qc-hint", icon("hand-pointer"), "Click a row to open the cell in the Inspector")),
    div(class = "qc-card", DT::DTOutput("tbl")))

  inspector <- shiny::tagList(
    div(class = "qc-toolbar",
        div(class = "qc-select-cell",
            shiny::selectizeInput("pick", NULL, choices = NULL, width = "260px",
                                  options = list(placeholder = paste0("Search cell", .qc_s[["ell"]])))),
        shiny::actionButton("prev", NULL, icon = icon("chevron-left"), class = "qc-btn qc-btn-icon",
                            title = paste0("Previous cell (", .qc_s[["larr"]], ")")),
        shiny::actionButton("nxt", NULL, icon = icon("chevron-right"), class = "qc-btn qc-btn-icon",
                            title = paste0("Next cell (", .qc_s[["rarr"]], ")")),
        shiny::uiOutput("insp_pos", inline = TRUE, class = "qc-pos"),
        div(class = "qc-spacer"),
        field("Browse", .qc_segmented("browse", c("All" = "all", "Failing" = "fail", "Passing" = "pass"),
                                      "all", "Cells to browse")),
        field("Decision", .qc_segmented("ovr", c("Auto" = "auto", "Keep" = "keep", "Exclude" = "exclude"),
                                        "auto", "Manual decision for this cell"))),
    div(class = "qc-card",
        shiny::uiOutput("insp_head", class = "qc-insp-head"),
        div(class = "qc-cn-row", .qc_cn_legend()),
        plotly::plotlyOutput("insp_profile", height = "340px")),
    div(class = "qc-card",
        card_head("Metrics", paste0("where this cell sits among all cells ", .qc_s[["mid"]], " the line marks this cell")),
        shiny::uiOutput("insp_tiles")))

  main <- tags$main(
    class = "qc-main",
    shiny::uiOutput("tab_status", class = "qc-tab-status"),
    bslib::navset_pill(
      id = "tab",
      bslib::nav_panel(tab_title("gauge-high", "Overview"), value = "overview", overview),
      bslib::nav_panel(tab_title("braille", "Explore"), value = "explore", explore),
      bslib::nav_panel(tab_title("table-list", "Cells"), value = "cells", cells),
      bslib::nav_panel(tab_title("dna", "Inspector"), value = "inspector", inspector)))

  bslib::page(
    title = paste0("ASCAT.sc ", .qc_s[["mid"]], " Quality control"), lang = "en", theme = .qc_bs_theme(),
    # Theme before first paint, so dark mode never flashes white.
    tags$head(tags$script(shiny::HTML(sprintf(
      "document.documentElement.setAttribute('data-bs-theme', '%s');", theme))),
      tags$link(rel = "icon", type = "image/webp", href = .qc_logo_uri("dark"))),
    # In <body>, after every stylesheet Shiny adds to <head>, so it wins ties.
    tags$style(shiny::HTML(.qc_css())),
    tags$script(shiny::HTML(.qc_js)),
    div(class = "qc-app", header, div(class = "qc-body", sidebar, main)))
}

# ── Server ───────────────────────────────────────────────────────────────────
.qc_server <- function(qd, res, st, theme, label = "res")
{
  M  <- qd$metrics
  ms <- M$metric
  fl <- qd$flags$flag
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  tags <- shiny::tags; div <- shiny::div; span <- shiny::span; icon <- shiny::icon

  # Profiles are extracted on demand and cached per cell. If extraction fails
  # (e.g. an unusual res layout) the reason is shown instead of the plot.
  prof_cache <- new.env(parent = emptyenv())
  prof_err <- new.env(parent = emptyenv())
  get_profile <- function(id) {
    if (!qd$has_profiles) return(NULL)
    if (exists(id, envir = prof_cache, inherits = FALSE)) return(get(id, envir = prof_cache))
    j <- qd$prof_idx[match(id, qd$ids)]
    sex <- if (!is.na(j) && length(res$sex) >= j) res$sex[j] else NULL
    pd <- tryCatch(.qc_profile_data(res, j, sex), error = function(e) {
      assign(id, conditionMessage(e), envir = prof_err); NULL })
    if (length(ls(prof_cache)) > 80) rm(list = ls(prof_cache), envir = prof_cache)
    assign(id, pd, envir = prof_cache)
    pd
  }
  no_profile_msg <- function(id) {
    err <- get0(id, envir = prof_err, inherits = FALSE)
    if (is.null(err)) "No copy-number profile is available for this cell."
    else paste("Could not build this cell's profile:", err)
  }
  status_pill <- function(s) {
    lab <- c(pass = "Pass", fail = "Fail", kept = paste0("Kept ", .qc_s[["mid"]], " manual"),
             excluded = paste0("Excluded ", .qc_s[["mid"]], " manual"))[[s]]
    ic <- if (s %in% c("pass", "kept")) "circle-check" else "circle-xmark"
    span(class = paste("qc-status", s), icon(ic), lab)
  }
  unescape <- function(x) {
    for (p in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&#39;", "'"), c("&amp;", "&")))
      x <- gsub(p[1], p[2], x, fixed = TRUE)
    x
  }

  ev <- function(event, source) plotly::event_data(event, source = source)

  function(input, output, session)
  {
    # The Explore plot is not drawn until its tab is opened, and plotly warns
    # (after each flush) about events of plots it has not seen yet. Register
    # them up front; the IDs follow plotly's "<event>-<source>" convention.
    session$userData$plotlyShinyEventIDs <- unique(c(
      session$userData$plotlyShinyEventIDs,
      "plotly_click-explore", "plotly_selected-explore", "plotly_deselect-explore",
      "plotly_click-corr"))

    pal <- shiny::reactive(.qc_plot_palette(
      if (input$dark_mode %in% c("dark", "light")) input$dark_mode else theme))

    # ---- Gate state --------------------------------------------------------
    gate_b <- function(m) {
      v <- input[[paste0("gate_", m)]]
      if (length(v) == 2) v else st$state$bounds[[m]]
    }
    gate_on <- function(m) {
      v <- input[[paste0("on_", m)]]
      if (is.null(v)) isTRUE(st$state$on[[m]]) else isTRUE(v)
    }
    bounds   <- shiny::reactive(stats::setNames(lapply(ms, gate_b), ms))
    on       <- shiny::reactive(stats::setNames(vapply(ms, gate_on, logical(1)), ms))
    flags_on <- shiny::reactive(stats::setNames(vapply(fl, function(f) {
      v <- input[[paste0("flag_", f)]]
      if (is.null(v)) isTRUE(st$state$flags_on[[f]]) else isTRUE(v)
    }, logical(1)), fl))
    overrides <- shiny::reactiveVal(st$state$overrides)
    k_val <- shiny::reactive({
      k <- suppressWarnings(as.numeric(input$k))
      if (length(k) != 1 || !is.finite(k) || k <= 0) st$state$k else k
    })

    cls <- shiny::reactive(.qc_classify(qd, bounds(), on(), flags_on(), overrides()))
    # Plots redraw once dragging pauses
    gates_d <- shiny::debounce(shiny::reactive(
      list(cls = cls(), bounds = bounds(), on = on())), 150)

    state_now <- shiny::reactive(list(bounds = bounds(), on = on(), flags_on = flags_on(),
                                      overrides = overrides(), k = k_val()))
    shiny::observe(st$state <- state_now())

    set_override <- function(ids, value) {
      ov <- overrides()
      if (value == "auto") ov <- ov[!names(ov) %in% ids] else ov[ids] <- value
      overrides(ov)
    }

    # ---- Sidebar gates -----------------------------------------------------
    lapply(ms, function(m) {
      sp <- M[m, ]
      med <- .qc_fmt(if (sp$log) 10^sp$median else sp$median)
      output[[paste0("gs_", m)]] <- shiny::renderUI({
        b <- gate_b(m); is_on <- gate_on(m)
        nfail <- if (is_on) sum(.qc_gate_fail(qd$S[, m], sp, b)) else 0L
        shiny::tagList(
          span(class = paste("qc-gate-badge", if (nfail == 0) "is-zero"),
               title = "Cells outside this threshold", paste(fmt(nfail), "fail")),
          div(class = "qc-gate-read",
              span(class = "qc-gate-limit", .qc_gate_text(sp, b)),
              span(class = "qc-gate-med", paste("median", med))),
          .qc_spark_svg(qd$hist[[m]], b[1], b[2], on = is_on))
      })
    })

    shiny::observeEvent(input$apply_k, {
      k <- k_val()
      for (m in ms) if (gate_on(m))
        shiny::updateSliderInput(session, paste0("gate_", m), value = .qc_default_bounds(M[m, ], k))
      shiny::showNotification(sprintf(paste0("Thresholds in use set to median ", .qc_s[["pm"]], " %s ", .qc_s[["times"]], " MAD."), k), duration = 3)
    })
    shiny::observeEvent(input$reset, {
      d <- st$defaults
      for (m in ms) {
        shiny::updateSliderInput(session, paste0("gate_", m), value = d$bounds[[m]])
        bslib::update_switch(paste0("on_", m), value = isTRUE(d$on[[m]]))
      }
      for (f in fl) bslib::update_switch(paste0("flag_", f), value = FALSE)
      shiny::updateNumericInput(session, "k", value = d$k)
      shiny::showNotification("Thresholds reset to their defaults. Manual decisions were kept.", duration = 3)
    })

    # ---- Summary ------------------------------------------------------------
    output$tab_status <- shiny::renderUI({
      np <- sum(cls()$pass)
      shiny::tagList(span(class = "qc-dot pass"),
                     span(tags$b(fmt(np)), " of ", fmt(qd$n), " cells pass"),
                     span(class = "qc-faint", .qc_pct(np, qd$n)))
    })

    output$kpis <- shiny::renderUI({
      cl <- cls(); n <- qd$n; np <- sum(cl$pass); nf <- n - np
      ov <- overrides(); nk <- sum(ov == "keep"); ne <- sum(ov == "exclude")
      kpi <- function(ic, label, value, small = NULL, sub = NULL, meter = NULL)
        div(class = "qc-kpi",
            div(class = "qc-kpi-label", ic, label),
            div(class = "qc-kpi-value", value, if (!is.null(small)) tags$small(small)),
            if (!is.null(sub)) div(class = "qc-kpi-sub", sub),
            if (!is.null(meter)) div(class = "qc-meter", div(style = sprintf("width:%.1f%%", 100 * meter))))
      manual <- if (nk + ne == 0) "No manual decisions" else shiny::tagList(
        sprintf(paste0("Manual: %d kept ", .qc_s[["mid"]], " %d excluded ", .qc_s[["mid"]], " "), nk, ne),
        shiny::actionLink("clear_ovr", "clear"))
      div(class = "qc-kpis",
          kpi(icon("layer-group"), "Cells", fmt(n),
              sub = if (is.finite(qd$n_bins)) paste(fmt(qd$n_bins), "bins per cell") else "in res$QC_metrics"),
          kpi(icon("circle-check", class = "icon-pass"), "Passing", fmt(np), .qc_pct(np, n),
              meter = if (n) np / n else 0),
          kpi(icon("circle-xmark", class = "icon-fail"), "Failing", fmt(nf), .qc_pct(nf, n),
              sub = sprintf("%s fail a metric or flag", fmt(sum(!cl$gate_pass)))),
          kpi(icon("sliders"), "Metrics in use", sum(on()), sprintf("of %d", length(ms)),
              sub = manual))
    })
    shiny::observeEvent(input$clear_ovr, {
      ov <- overrides()
      shiny::showModal(shiny::modalDialog(
        title = "Clear manual decisions?",
        sprintf("This removes %d keep and %d exclude decisions. Thresholds are not changed.",
                sum(ov == "keep"), sum(ov == "exclude")),
        easyClose = TRUE,
        footer = shiny::tagList(shiny::modalButton("Cancel"),
                                shiny::actionButton("clear_ovr_ok", "Clear", class = "qc-btn qc-btn-primary"))))
    })
    shiny::observeEvent(input$clear_ovr_ok, { overrides(stats::setNames(character(), character())); shiny::removeModal() })

    # Pass / fail per metric (and flag) in use; the last row combines them all.
    output$impact <- shiny::renderUI({
      cl <- cls(); o <- on(); fo <- flags_on()
      F <- cbind(cl$fails[, o, drop = FALSE], cl$flag_fails[, fo, drop = FALSE])
      if (ncol(F) == 0)
        return(div(class = "qc-empty", icon("sliders"), "No metric is in use. Switch one on in the panel."))
      n <- qd$n
      labs <- c(M$label[o], qd$flags$label[fo])
      n_fail <- colSums(F)
      pct <- function(x) .qc_pct(x, n)
      row <- function(label, nf, total = FALSE) {
        np <- n - nf
        tip <- sprintf("%s: %s pass (%s), %s fail (%s)", label, fmt(np), pct(np), fmt(nf), pct(nf))
        extra <- if (total) "is-total"
        shiny::tagList(
          div(class = paste("qc-pf-label", extra), title = tip, label),
          div(class = paste("qc-pf-bar", extra), title = tip,
              if (np > 0) div(class = "qc-pf-pass", style = sprintf("flex: %d 1 0", np)),
              if (nf > 0) div(class = "qc-pf-fail", style = sprintf("flex: %d 1 0", nf))),
          div(class = paste("qc-pf-n", extra), title = tip, fmt(np)),
          div(class = paste("qc-pf-n", extra, if (nf == 0) "is-zero"), title = tip, fmt(nf)))
      }
      shiny::tagList(
        div(class = "qc-pf",
            div(class = "qc-pf-h", "Metric"),
            div(class = "qc-pf-h"),
            div(class = "qc-pf-h qc-pf-hn", span(class = "qc-dot pass"), "Passing"),
            div(class = "qc-pf-h qc-pf-hn", span(class = "qc-dot fail"), "Failing"),
            lapply(order(-n_fail), function(j) row(labs[j], n_fail[j])),
            div(class = "qc-pf-sep"),
            row("All metrics", sum(rowSums(F) > 0), total = TRUE)),
        div(class = "qc-hint", style = "margin-top:12px",
            "All metrics: a cell passes only if it is within every threshold.",
            if (length(overrides())) " Manual keep / exclude decisions are not counted here."))
    })

    output$corr <- plotly::renderPlotly(.qc_plot_corr(qd, pal()))
    shiny::observeEvent(ev("plotly_click", "corr"), {
      d <- ev("plotly_click", "corr")
      xm <- M$metric[match(d$x, M$label)]; ym <- M$metric[match(d$y, M$label)]
      if (length(xm) != 1 || is.na(xm) || is.na(ym) || xm == ym) return()
      shiny::updateSelectInput(session, "xvar", selected = xm)
      shiny::updateSelectInput(session, "yvar", selected = ym)
      bslib::nav_select("tab", "explore")
    })

    lapply(ms, function(m) {
      output[[paste0("hist_", m)]] <- plotly::renderPlotly({
        g <- gates_d()
        .qc_plot_hist(qd, m, g$bounds[[m]], g$on[[m]], g$cls$pass, pal())
      })
      output[[paste0("distgate_", m)]] <- shiny::renderText(
        if (gate_on(m)) .qc_gate_text(M[m, ], gate_b(m)) else "not used")
    })

    # ---- Explore ------------------------------------------------------------
    sel <- shiny::reactiveVal(NULL)
    lasso <- shiny::reactiveVal(character())

    output$explore <- plotly::renderPlotly({
      g <- gates_d()
      shiny::req(input$xvar, input$yvar)
      dim <- if (identical(input$dim, "2d")) "2d" else "3d"
      if (dim == "3d") shiny::req(input$zvar)
      .qc_plot_scatter(qd, input$xvar, input$yvar, input$zvar, g$cls$status, pal(),
                       g$bounds, g$on, selected = sel(),
                       color_by = if (isTRUE(input$color_by %in% c("status", ms))) input$color_by else "status",
                       dim = dim)
    })
    shiny::observeEvent(ev("plotly_click", "explore"), {
      id <- ev("plotly_click", "explore")$customdata
      id <- if (length(id)) as.character(unlist(id))[1] else NULL
      if (!is.null(id) && id %in% qd$ids) sel(id)
    })
    shiny::observeEvent(ev("plotly_selected", "explore"), {
      d <- ev("plotly_selected", "explore")
      ids <- if (is.null(d) || !"customdata" %in% names(d)) character()
             else intersect(unique(as.character(unlist(d$customdata))), qd$ids)
      lasso(ids)
    }, ignoreNULL = FALSE)
    shiny::observeEvent(ev("plotly_deselect", "explore"), lasso(character()))

    output$selbar <- shiny::renderUI({
      ids <- lasso()
      if (!length(ids)) return(NULL)
      s <- cls()$status[ids]
      div(class = "qc-selbar",
          icon("object-group"), span(tags$b(fmt(length(ids))), " cells selected"),
          span(class = "qc-muted", sprintf(paste0("%s pass ", .qc_s[["mid"]], " %s fail"), fmt(sum(s %in% c("pass", "kept"))),
                                           fmt(sum(s %in% c("fail", "excluded"))))),
          div(class = "qc-spacer"),
          shiny::actionButton("sel_keep", "Keep", icon = icon("circle-check"), class = "qc-btn qc-btn-sm"),
          shiny::actionButton("sel_excl", "Exclude", icon = icon("circle-xmark"), class = "qc-btn qc-btn-sm qc-btn-danger"),
          shiny::actionButton("sel_auto", "Back to auto", icon = icon("rotate-left"), class = "qc-btn qc-btn-sm"),
          shiny::actionButton("sel_clear", NULL, icon = icon("xmark"), class = "qc-btn qc-btn-sm qc-btn-ghost",
                              title = "Clear the selection"))
    })
    shiny::observeEvent(input$sel_keep, set_override(lasso(), "keep"))
    shiny::observeEvent(input$sel_excl, set_override(lasso(), "exclude"))
    shiny::observeEvent(input$sel_auto, set_override(lasso(), "auto"))
    shiny::observeEvent(input$sel_clear, lasso(character()))

    output$has_sel <- shiny::reactive(!is.null(sel()))
    shiny::outputOptions(output, "has_sel", suspendWhenHidden = FALSE)

    cell_summary <- function(id) {
      i <- match(id, qd$ids); cl <- cls()
      fails <- .qc_failed_labels(list(fails = cl$fails[i, , drop = FALSE],
                                      flag_fails = cl$flag_fails[i, , drop = FALSE]), list(
                                        metrics = M, flags = qd$flags, n = 1))
      meta <- list()
      if (!is.null(qd$df$ploidy)) meta <- c(meta, list("ploidy ", tags$b(.qc_fmt(qd$df$ploidy[i]))))
      if (!is.null(qd$df$purity)) meta <- c(meta, list(paste0(" ", .qc_s[["mid"]], " purity "), tags$b(.qc_fmt(qd$df$purity[i]))))
      if (!is.null(qd$df$n_reads_total))
        meta <- c(meta, list(paste0(" ", .qc_s[["mid"]], " "), tags$b(.qc_fmt(qd$df$n_reads_total[i])), " reads"))
      list(id = id, status = cl$status[[i]], fails = fails, meta = meta)
    }

    output$explore_sel_head <- shiny::renderUI({
      shiny::req(sel())
      s <- cell_summary(sel())
      shiny::tagList(span(class = "qc-cell-id", s$id), status_pill(s$status),
                     span(class = "qc-cell-meta", s$meta),
                     if (nzchar(s$fails)) span(class = "qc-cell-fails", "Fails: ", tags$b(s$fails)))
    })
    output$explore_profile <- plotly::renderPlotly({
      shiny::req(sel())
      shiny::validate(shiny::need(qd$has_profiles,
        "Copy-number profiles need res$allTracks.processed, res$allProfiles and res$allSolutions."))
      pd <- get_profile(sel())
      shiny::validate(shiny::need(!is.null(pd), no_profile_msg(sel())))
      .qc_plot_profile(pd, pal(), compact = TRUE)
    })
    shiny::observeEvent(input$to_inspector, bslib::nav_select("tab", "inspector"))

    # ---- Cells table -------------------------------------------------------
    tbl_data <- shiny::reactive({
      cl <- cls(); ov <- overrides()
      keep <- switch(if (is.null(input$tbl_filter)) "all" else input$tbl_filter,
                     pass = cl$pass, fail = !cl$pass, manual = qd$ids %in% names(ov),
                     rep(TRUE, qd$n))
      lab <- c(pass = "Pass", fail = "Fail", kept = "Kept", excluded = "Excluded")
      # First two failed gates, "+N" for the rest, full list on hover
      fails <- strsplit(.qc_failed_labels(cl, qd), ", ", fixed = TRUE)
      fails_html <- vapply(fails, function(f) {
        if (!length(f)) return("")
        txt <- htmltools::htmlEscape(paste(utils::head(f, 2), collapse = ", "))
        if (length(f) > 2) txt <- sprintf('%s <span class="qc-more">+%d</span>', txt, length(f) - 2)
        sprintf('<span class="qc-fails" title="%s">%s</span>',
                htmltools::htmlEscape(paste(f, collapse = ", "), attribute = TRUE), txt)
      }, character(1))
      d <- data.frame(Cell = qd$ids,
                      Status = sprintf('<span class="qc-status %s">%s</span>', cl$status, lab[cl$status]),
                      `Failed metrics` = fails_html,
                      check.names = FALSE, stringsAsFactors = FALSE)
      for (m in ms) d[[M[m, "label"]]] <- qd$df[[m]]
      d[keep, , drop = FALSE]
    })
    output$tbl <- DT::renderDT({
      d <- shiny::isolate(tbl_data())
      dt <- DT::datatable(
        d, rownames = FALSE, escape = c(-2, -3), selection = "none", style = "bootstrap5",
        class = "table table-sm", width = "100%",
        options = list(pageLength = 25, lengthMenu = c(25, 50, 100, 500), scrollX = TRUE,
                       dom = "<'d-flex justify-content-between align-items-center mb-2'fl>t<'d-flex justify-content-between align-items-center mt-2'ip>",
                       order = list(), autoWidth = FALSE,
                       language = list(search = "", searchPlaceholder = paste0("Search cells", .qc_s[["ell"]]),
                                       lengthMenu = "_MENU_ per page",
                                       info = paste0("_START_", .qc_s[["ndash"]], "_END_ of _TOTAL_"),
                                       infoFiltered = "(of _MAX_)", zeroRecords = "No cells match")),
        callback = DT::JS(
          "table.on('click', 'tbody tr', function () {",
          "  var d = table.row(this).data();",
          "  if (d) Shiny.setInputValue('tbl_click', d[0], {priority: 'event'});",
          "});"))
      lab_int <- M$label[M$integer]
      lab_big <- M$label[M$log]
      lab_dec <- setdiff(M$label, c(lab_int, lab_big))
      if (length(lab_dec)) dt <- DT::formatSignif(dt, lab_dec, digits = 3)
      if (length(c(lab_int, lab_big))) dt <- DT::formatRound(dt, c(lab_int, lab_big), digits = 0)
      dt
    }, server = TRUE)
    tbl_proxy <- DT::dataTableProxy("tbl")
    shiny::observeEvent(tbl_data(), DT::replaceData(tbl_proxy, tbl_data(), resetPaging = FALSE,
                                                    rownames = FALSE), ignoreInit = TRUE)
    shiny::observeEvent(input$tbl_click, {
      id <- input$tbl_click
      if (!id %in% qd$ids) id <- unescape(id)
      if (id %in% qd$ids) { sel(id); bslib::nav_select("tab", "inspector") }
    })

    # ---- Inspector ---------------------------------------------------------
    browse_ids <- shiny::reactive({
      p <- cls()$pass
      switch(if (is.null(input$browse)) "all" else input$browse,
             fail = qd$ids[!p], pass = qd$ids[p], qd$ids)
    })
    insp_id <- shiny::reactive({
      id <- sel()
      if (!is.null(id)) return(id)
      b <- browse_ids()
      if (length(b)) b[1] else qd$ids[1]
    })
    go <- function(dir) {
      b <- browse_ids()
      if (!length(b)) return()
      i <- match(insp_id(), b)
      j <- if (is.na(i)) (if (dir > 0) 1 else length(b)) else max(1, min(length(b), i + dir))
      sel(b[j])
    }
    shiny::observeEvent(input$prev, go(-1))
    shiny::observeEvent(input$nxt, go(1))
    shiny::observeEvent(input$qc_nav, if (input$tab %in% c("inspector", "explore")) go(input$qc_nav))
    shiny::observeEvent(input$browse, {
      b <- browse_ids()
      if (length(b) && (is.null(sel()) || !sel() %in% b)) sel(b[1])
    }, ignoreInit = TRUE)
    # The cell shown in the Inspector is the selected cell everywhere else too
    shiny::observeEvent(input$tab, {
      if (identical(input$tab, "inspector") && is.null(sel())) sel(insp_id())
    })

    shiny::updateSelectizeInput(session, "pick", choices = qd$ids, server = TRUE,
                                selected = shiny::isolate(insp_id()))
    shiny::observeEvent(input$pick, {
      if (length(input$pick) && input$pick %in% qd$ids && !identical(input$pick, insp_id()))
        sel(input$pick)
    }, ignoreInit = TRUE)
    shiny::observeEvent(insp_id(), {
      if (!identical(input$pick, insp_id()))
        shiny::updateSelectizeInput(session, "pick", choices = qd$ids, server = TRUE,
                                    selected = insp_id())
    }, ignoreInit = TRUE)

    # Decision control mirrors the inspected cell's override
    shiny::observeEvent(list(insp_id(), overrides()), {
      cur <- overrides()[insp_id()]
      cur <- if (is.na(cur)) "auto" else unname(cur)
      if (!identical(input$ovr, cur)) shiny::updateRadioButtons(session, "ovr", selected = cur)
    })
    shiny::observeEvent(input$ovr, {
      cur <- overrides()[insp_id()]
      cur <- if (is.na(cur)) "auto" else unname(cur)
      if (!identical(cur, input$ovr)) set_override(insp_id(), input$ovr)
    }, ignoreInit = TRUE)

    output$insp_pos <- shiny::renderUI({
      b <- browse_ids(); i <- match(insp_id(), b)
      what <- switch(if (is.null(input$browse)) "all" else input$browse,
                     fail = "failing", pass = "passing", "cells")
      if (!length(b)) return(span(sprintf("no %s cells", what)))
      span(if (is.na(i)) sprintf("not in the %s list (%s)", what, fmt(length(b)))
           else sprintf("%s of %s %s", fmt(i), fmt(length(b)), what))
    })

    output$insp_head <- shiny::renderUI({
      s <- cell_summary(insp_id())
      pd <- get_profile(insp_id())
      shiny::tagList(
        div(class = "qc-cell-head",
            span(class = "qc-cell-id", s$id), status_pill(s$status),
            span(class = "qc-cell-meta", s$meta),
            if (!is.null(pd)) span(class = "qc-cell-meta qc-faint",
                                   sprintf(paste0(.qc_s[["mid"]], " %s segments"), fmt(nrow(pd$segs))))),
        div(class = "qc-cell-fails",
            if (nzchar(s$fails)) shiny::tagList("Fails: ", tags$b(s$fails))
            else if (s$status == "excluded") "Passes every metric, excluded manually."
            else "Passes every metric in use."))
    })
    output$insp_profile <- plotly::renderPlotly({
      shiny::validate(shiny::need(qd$has_profiles,
        "Copy-number profiles need res$allTracks.processed, res$allProfiles and res$allSolutions."))
      pd <- get_profile(insp_id())
      shiny::validate(shiny::need(!is.null(pd), no_profile_msg(insp_id())))
      .qc_plot_profile(pd, pal())
    })

    output$insp_tiles <- shiny::renderUI({
      i <- match(insp_id(), qd$ids); cl <- cls(); b <- bounds(); o <- on()
      tiles <- lapply(ms, function(m) {
        sp <- M[m, ]; s <- qd$S[i, m]
        pct <- if (is.finite(s)) round(100 * mean(qd$S[, m] <= s, na.rm = TRUE)) else NA
        failing <- isTRUE(cl$fails[i, m])
        st_ <- if (!o[[m]]) c("", "off", "not used")
               else if (failing) c("fail", "fail", paste("outside", .qc_gate_text(sp, b[[m]])))
               else c("pass", "pass", paste("within", .qc_gate_text(sp, b[[m]])))
        div(class = paste("qc-tile", if (o[[m]] && failing) "is-fail"),
            div(class = "qc-tile-label", title = sp$label, sp$label),
            div(class = "qc-tile-row",
                span(class = "qc-tile-value", .qc_fmt(qd$df[[m]][i])),
                if (!is.na(pct)) span(class = "qc-tile-pct", title = sprintf("Higher than or equal to %d%% of cells", pct),
                                      sprintf("P%d", pct))),
            .qc_spark_svg(qd$hist[[m]], b[[m]][1], b[[m]][2], on = o[[m]], marker = s, height = 22),
            div(class = paste("qc-tile-status", st_[1]), span(class = paste("qc-dot", st_[2])), st_[3]))
      })
      flag_tiles <- lapply(fl, function(f) {
        fr <- qd$flags[f, ]; v <- qd$df[[f]][i]
        hit <- isTRUE(v %in% fr$exclude_if); act <- isTRUE(flags_on()[[f]])
        st_ <- if (!act) c("", "off", "flag off") else if (hit) c("fail", "fail", "excluded by flag")
               else c("pass", "pass", "not flagged")
        div(class = paste("qc-tile", if (act && hit) "is-fail"),
            div(class = "qc-tile-label", title = fr$label, fr$label),
            div(class = "qc-tile-row", span(class = "qc-tile-value",
                                            if (is.na(v)) .qc_s[["dash"]] else if (isTRUE(v)) "Yes" else "No")),
            div(class = paste("qc-tile-status", st_[1]), style = "margin-top:6px",
                span(class = paste("qc-dot", st_[2])), st_[3]))
      })
      div(class = "qc-tiles", tiles, flag_tiles)
    })

    # ---- Export & finish ---------------------------------------------------
    result_now <- function() .qc_result(qd, state_now())
    stamp <- function() format(Sys.Date())
    output$dl_pass <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_passing_%s.csv", stamp()),
      content = function(file) {
        t <- result_now()$cells
        utils::write.csv(t[t$qc_pass, , drop = FALSE], file, row.names = FALSE)
      })
    output$dl_fail <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_failing_%s.csv", stamp()),
      content = function(file) {
        t <- result_now()$cells
        utils::write.csv(t[!t$qc_pass, , drop = FALSE], file, row.names = FALSE)
      })
    output$dl_table <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_all_cells_%s.csv", stamp()),
      content = function(file) utils::write.csv(result_now()$cells, file, row.names = FALSE))
    output$dl_gates <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_thresholds_%s.csv", stamp()),
      content = function(file) {
        r <- result_now()
        g <- data.frame(type = "metric", r$thresholds, exclude_if = NA, mad_k = r$k)
        f <- r$flags
        if (nrow(f)) {
          f <- data.frame(type = "flag", metric = f$flag, label = f$label, enabled = f$enabled,
                          direction = NA, min = NA, max = NA, data_min = NA, data_max = NA,
                          n_fail = f$n_fail, exclude_if = f$exclude_if, mad_k = r$k)
          g <- rbind(g, f)
        }
        utils::write.csv(g, file, row.names = FALSE)
      })

    shiny::observeEvent(input$done, {
      np <- sum(cls()$pass)
      shiny::insertUI("body", "beforeEnd", immediate = TRUE, div(
        class = "qc-done", role = "status",
        div(class = "qc-done-card", .qc_logo(),
            div(class = "qc-done-title", "Quality control finished"),
            div(class = "qc-done-text",
                sprintf("%s of %s cells pass (%s). The result is saved in %s in your R session; you can close this tab.",
                        fmt(np), fmt(qd$n), .qc_pct(np, qd$n), label)),
            tags$code(paste(paste0(label, c("$pass", "$fail", "$cells")),
                            collapse = paste0(" ", .qc_s[["mid"]], " "))))))
      shiny::stopApp()
    })
  }
}

# =============================================================================
# Internal helpers. Names start with "." so that exportPattern("^[[:alpha:]]+")
# in NAMESPACE keeps them internal (and hidden from ls() when sourced).
# =============================================================================

# ── Metric & flag metadata ───────────────────────────────────────────────────
# dir:  "upper" = high values are bad (gate has a max), "lower" = low values
#       are bad (gate has a min), "both" = two-sided.
# log:    gate and plot on a log10 scale.
# on:     gate enabled by default.
# robust: slider covers median +/- 20 MAD rather than the full data range, so
#         a few extreme cells do not squeeze the useful thresholds together.
.qc_metric_meta <- function()
{
  m <- function(label, group, dir, help, on = TRUE, log = FALSE, robust = group %in% c("Noise", "Fit"))
    list(label = label, group = group, dir = dir, help = help, on = on, log = log, robust = robust)
  list(
    n_reads_total = m("Total reads", "Coverage", "lower", log = TRUE,
      "Total number of reads across all bins. Low counts give noisy copy-number profiles."),
    MAPD_gc_corrected = m("MAPD", "Noise", "upper",
      "Median absolute pairwise difference of GC-corrected logR between neighbouring bins. Measures bin-level noise."),
    spikiness = m("Spikiness", "Noise", "upper",
      "Median second-order difference of logR. Detects isolated outlier bins. Lower is smoother."),
    MBRSM_dispersion = m("Segment dispersion", "Noise", "upper",
      "Scaled MAD of within-segment residuals. Measures intra-segment noise."),
    autocorr = m("Autocorrelation", "Noise", "upper",
      "Lag-1 autocorrelation of segment residuals. High values suggest wave-like bias (GC, replication timing)."),
    sd_noise = m("Noise (getFilters)", "Noise", "upper",
      "Median absolute difference of consecutive logR bins, as computed by getFilters()."),
    distance_integer_logR = m("LogR integer distance", "Fit", "upper",
      "Weighted distance of segment copy numbers from the nearest integer. Lower is a better fit."),
    distance_integer_BAF = m("BAF integer distance", "Fit", "upper",
      "Weighted distance of segment BAFs from the closest integer allelic state. Lower is cleaner."),
    MSRSI_non_integerness = m("Non-integerness", "Fit", "upper",
      "Scaled MAD of segment logR deviations from the integer-expected logR (MSRSI). Lower is more integer-like."),
    breakpoints = m("Breakpoints", "Fit", "upper",
      "Number of copy-number state changes across the genome. Very high counts usually mean noise is being segmented."),
    state_mode = m("Modal CN state", "Solution", "both", on = FALSE,
      "Most frequent copy-number state. Typically 2 for near-diploid cells."),
    ploidy = m("Ploidy", "Solution", "both", on = FALSE,
      "Fitted ploidy from the purity/ploidy grid search (res$allSolutions)."),
    purity = m("Purity", "Solution", "both", on = FALSE,
      "Fitted purity from the grid search (res$allSolutions). Single cells are expected close to 1.")
  )
}

.qc_group_order <- c("Coverage", "Noise", "Fit", "Solution", "Other")

# Logical columns become exclusion switches. exclude_if is the value that
# fails a cell when the switch is on.
.qc_flag_meta <- function()
{
  f <- function(label, exclude_if, help)
    list(label = label, exclude_if = exclude_if, help = help)
  list(
    ambiguous_flag = f("Ambiguous fit", TRUE,
      "searchGrid() could not settle on a purity/ploidy solution."),
    ambiguous = f("Ambiguous fit", TRUE,
      "searchGrid() could not settle on a purity/ploidy solution."),
    noise_filter_out = f("Noise outlier (getFilters)", TRUE,
      "Flagged by getFilters() as being among the noisiest cells."),
    n_reads_filter_out = f("Low reads (getFilters)", TRUE,
      "Flagged by getFilters() as being among the cells with fewest reads."),
    keep_cell_filter = f("Not kept by getFilters()", FALSE,
      "getFilters() did not keep this cell (keep_cell_filter is FALSE).")
  )
}

# ── Small utilities ──────────────────────────────────────────────────────────
.qc_or <- function(a, b) if (is.null(a)) b else a

# Typographic symbols used in labels. They are built from code points, not
# written into strings: R code must be ASCII to be portable, and in a C locale
# (common on clusters) R turns such characters into the text "<U+2264>".
.qc_s <- vapply(c(le = 0x2264, ge = 0x2265, dash = 0x2014, ndash = 0x2013,
                  pm = 0x00b1, times = 0x00d7, mid = 0x00b7, ell = 0x2026,
                  larr = 0x2190, rarr = 0x2192, uarr = 0x2191, darr = 0x2193,
                  lsaquo = 0x2039, rsaquo = 0x203a, rho = 0x03c1),
                intToUtf8, character(1))

# A "nice" slider step: 1, 2 or 5 times a power of ten.
.qc_nice_step <- function(x)
{
  if (!is.finite(x) || x <= 0) return(1)
  p <- 10^floor(log10(x))
  m <- x / p
  (if (m < 1.5) 1 else if (m < 3.5) 2 else if (m < 7.5) 5 else 10) * p
}

# Snap to the step grid; the rounding removes floating-point dust.
.qc_snap <- function(v, step, fn = round)
{
  digits <- max(0, -floor(log10(step))) + 2
  round(fn(v / step + 1e-9) * step, digits)
}

.qc_sig <- function(v, digits = 3)
{
  s <- formatC(signif(v, digits), digits = digits, format = "fg", flag = "#")
  sub("\\.$", "", trimws(s))
}

# Compact, human formatting: 1.29M, 65.3k, 0.243, 12. Vectorised, since hover
# labels format every cell.
.qc_fmt <- function(x, digits = 3)
{
  x <- as.numeric(x)
  out <- rep(.qc_s[["dash"]], length(x))
  ok <- is.finite(x)
  v <- x[ok]; a <- abs(v)
  mil <- a >= 1e6
  kil <- a >= 1e4 & !mil
  int <- !mil & !kil & v == round(v)
  oth <- !(mil | kil | int)
  r <- character(length(v))
  if (any(mil)) r[mil] <- paste0(.qc_sig(v[mil] / 1e6, digits), "M")
  if (any(kil)) r[kil] <- paste0(.qc_sig(v[kil] / 1e3, digits), "k")
  if (any(int)) r[int] <- format(v[int], big.mark = ",", scientific = FALSE, trim = TRUE)
  if (any(oth)) r[oth] <- .qc_sig(v[oth], digits)
  out[ok] <- r
  out
}

.qc_pct <- function(n, total)
  if (total == 0) "0%" else paste0(formatC(100 * n / total, format = "f", digits = 1), "%")

# Histogram counts over fixed breaks; values outside are clamped into the
# first / last bin.
.qc_hist_counts <- function(s, breaks)
{
  s <- s[is.finite(s)]
  s <- pmin(pmax(s, breaks[1]), breaks[length(breaks)])
  graphics::hist(s, breaks = breaks, plot = FALSE)$counts
}

# ── Data preparation ─────────────────────────────────────────────────────────
# Returns everything the app needs, computed once:
#   df      the QC table with a cell_id column (+ ploidy/purity if available)
#   ids     cell ids
#   S       numeric matrix of gate metrics on the slider scale (log10 if log)
#   metrics one row per gate metric: label, group, dir, log, slider range ...
#   flags   one row per logical flag column
#   hist    per-metric histogram (on the slider scale) for the sparklines
#   prof_idx index of each cell in res$allProfiles (NA if not found)
.qc_prepare <- function(res, previous = NULL)
{
  qc <- as.data.frame(res$QC_metrics, stringsAsFactors = FALSE)
  n  <- nrow(qc)
  if (n == 0) stop("res$QC_metrics has no rows.")

  rn <- rownames(qc)
  ids <- if (!is.null(rn) && !identical(rn, as.character(seq_len(n)))) rn
         else if (length(res$allProfiles) == n && !is.null(names(res$allProfiles))) names(res$allProfiles)
         else paste0("cell_", seq_len(n))
  ids <- make.unique(as.character(ids))

  # Map each QC row to its entry in res$allProfiles / allSolutions
  match_idx <- function(lst) {
    if (is.null(lst)) return(rep(NA_integer_, n))
    nm <- names(lst)
    if (!is.null(nm) && all(ids %in% nm)) return(match(ids, nm))
    if (length(lst) == n) return(seq_len(n))
    rep(NA_integer_, n)
  }
  prof_idx <- match_idx(res$allProfiles)
  sol_idx  <- match_idx(res$allSolutions)

  # Fitted ploidy / purity / ambiguity from the grid search, when available
  if (!all(is.na(sol_idx))) {
    sol_field <- function(f) vapply(sol_idx, function(j) {
      s <- if (is.na(j)) NULL else res$allSolutions[[j]]
      v <- if (is.list(s)) s[[f]] else NULL
      if (length(v) == 1 && !is.na(suppressWarnings(as.numeric(v)))) as.numeric(v) else NA_real_
    }, numeric(1))
    if (!"ploidy" %in% names(qc)) qc$ploidy <- sol_field("ploidy")
    if (!"purity" %in% names(qc)) qc$purity <- sol_field("purity")
    if (!any(c("ambiguous_flag", "ambiguous") %in% names(qc))) {
      amb <- sol_field("ambiguous")
      if (!all(is.na(amb))) qc$ambiguous <- as.logical(amb)
    }
  }

  meta_m <- .qc_metric_meta()
  meta_f <- .qc_flag_meta()

  # ---- gate metrics: numeric, not all NA, at least two distinct values ----
  is_num <- vapply(qc, function(x) is.numeric(x) && !is.logical(x), logical(1))
  cand <- names(qc)[is_num]
  cand <- cand[vapply(cand, function(col) {
    x <- qc[[col]]
    length(unique(x[is.finite(x)])) >= 2
  }, logical(1))]

  prev_g <- if (is.list(previous)) .qc_or(previous[["thresholds"]], previous[["gates"]]) else NULL
  if (!is.data.frame(prev_g)) prev_g <- NULL
  specs <- lapply(cand, function(col) {
    x  <- qc[[col]]
    mm <- .qc_or(meta_m[[col]],
                 list(label = gsub("_", " ", col), group = "Other", dir = "both",
                      help = NA_character_, on = TRUE, log = FALSE, robust = TRUE))
    ok  <- is.finite(x)
    log <- isTRUE(mm$log) && all(x[ok] > 0)
    s   <- if (log) log10(x) else x
    rng <- range(s[is.finite(s)])
    med <- stats::median(s, na.rm = TRUE)
    md  <- stats::mad(s, na.rm = TRUE)
    if (isTRUE(mm$robust) && is.finite(md) && md > 0) {
      if (mm$dir %in% c("lower", "both")) rng[1] <- max(rng[1], med - 20 * md)
      if (mm$dir %in% c("upper", "both")) rng[2] <- min(rng[2], med + 20 * md)
    }
    # Keep the finite limits of a previous run inside the slider range
    if (!is.null(prev_g) && col %in% prev_g$metric) {
      lim <- unlist(prev_g[prev_g$metric == col, c("min", "max")])
      lim <- lim[is.finite(lim) & (!log | lim > 0)]
      if (length(lim)) rng <- range(c(rng, if (log) log10(lim) else lim))
    }
    is_int <- !log && all(x[ok] == round(x[ok]))
    step <- .qc_nice_step(diff(rng) / 200)
    if (is_int) step <- max(1, round(step))
    lo <- .qc_snap(rng[1], step, floor)
    hi <- .qc_snap(rng[2], step, ceiling)
    if (hi <= lo) hi <- lo + step
    data.frame(metric = col, label = mm$label, group = mm$group, dir = mm$dir,
               help = .qc_or(mm$help, NA_character_), on = isTRUE(mm$on),
               log = log, integer = is_int, min = lo, max = hi, step = step,
               median = med, mad = md,
               n_below = sum(s < lo, na.rm = TRUE), n_above = sum(s > hi, na.rm = TRUE),
               n_na = sum(!ok), stringsAsFactors = FALSE)
  })
  metrics <- do.call(rbind, specs)
  if (is.null(metrics)) stop("No numeric QC metrics with variation were found in res$QC_metrics.")
  ord <- order(match(metrics$group, .qc_group_order), match(metrics$metric, names(meta_m)),
               na.last = TRUE)
  metrics <- metrics[ord, , drop = FALSE]
  rownames(metrics) <- metrics$metric

  S <- vapply(metrics$metric, function(col) {
    x <- qc[[col]]
    if (metrics[col, "log"]) suppressWarnings(log10(x)) else as.numeric(x)
  }, numeric(n))
  S <- matrix(S, nrow = n, dimnames = list(ids, metrics$metric))

  # ---- flags: logical columns ----
  lg <- names(qc)[vapply(qc, is.logical, logical(1))]
  lg <- lg[vapply(lg, function(col) any(!is.na(qc[[col]])), logical(1))]
  flags <- do.call(rbind, lapply(lg, function(col) {
    fm <- .qc_or(meta_f[[col]], list(label = gsub("_", " ", col), exclude_if = TRUE,
                                     help = NA_character_))
    data.frame(flag = col, label = fm$label, exclude_if = fm$exclude_if,
               help = .qc_or(fm$help, NA_character_),
               n = sum(qc[[col]] %in% fm$exclude_if), stringsAsFactors = FALSE)
  }))
  if (is.null(flags))
    flags <- data.frame(flag = character(), label = character(), exclude_if = logical(),
                        help = character(), n = integer(), stringsAsFactors = FALSE)
  rownames(flags) <- flags$flag

  # ---- sparkline histograms (fixed 36 bins over the slider range; cells
  #      beyond the range pile up in the end bins) ----
  hist <- lapply(stats::setNames(metrics$metric, metrics$metric), function(col) {
    sp <- metrics[col, ]
    br <- seq(sp$min, sp$max, length.out = 37)
    cnt <- .qc_hist_counts(S[, col], br)
    list(breaks = br, counts = cnt, mids = (br[-1] + br[-length(br)]) / 2)
  })

  df <- cbind(data.frame(cell_id = ids, stringsAsFactors = FALSE), qc)
  rownames(df) <- NULL

  has_profiles <- !is.null(res$allTracks.processed) && !is.null(res$allProfiles) &&
    !is.null(res$allSolutions) && !all(is.na(prof_idx))

  list(df = df, ids = ids, n = n, S = S, metrics = metrics, flags = flags,
       hist = hist, prof_idx = prof_idx, has_profiles = has_profiles,
       n_bins = if (has_profiles) tryCatch(sum(vapply(
         res$allTracks.processed[[prof_idx[which(!is.na(prof_idx))[1]]]]$lCTS,
         nrow, integer(1))), error = function(e) NA_integer_) else NA_integer_)
}

# Default gate on the slider scale: median +/- k * MAD, applied only on the
# side(s) where the metric's direction says values are bad.
.qc_default_bounds <- function(sp, k = 3)
{
  full <- c(sp$min, sp$max)
  if (!is.finite(sp$mad) || sp$mad == 0) return(full)
  lo <- if (sp$dir %in% c("lower", "both"))
    max(sp$min, .qc_snap(sp$median - k * sp$mad, sp$step, floor)) else sp$min
  hi <- if (sp$dir %in% c("upper", "both"))
    min(sp$max, .qc_snap(sp$median + k * sp$mad, sp$step, ceiling)) else sp$max
  c(lo, hi)
}

# Initial gate state: defaults, optionally overridden by a previous result.
# Gates that are off by default (ploidy, purity, ...) start fully open:
# those metrics are often multimodal, so median +/- MAD means little.
.qc_initial_state <- function(qd, k = 3, previous = NULL)
{
  M <- qd$metrics
  bounds <- lapply(stats::setNames(M$metric, M$metric), function(m)
    if (M[m, "on"]) .qc_default_bounds(M[m, ], k) else c(M[m, "min"], M[m, "max"]))
  on <- stats::setNames(M$on, M$metric)
  flags_on <- stats::setNames(rep(FALSE, nrow(qd$flags)), qd$flags$flag)
  overrides <- stats::setNames(character(), character())

  if (!is.null(previous)) {
    g <- .qc_or(previous[["thresholds"]], previous[["gates"]])   # $gates: older sessions
    if (is.data.frame(g) && all(c("metric", "enabled", "min", "max") %in% names(g))) {
      for (r in seq_len(nrow(g))) {
        m <- g$metric[r]
        if (!m %in% M$metric) next
        to_s <- function(v) if (M[m, "log"]) log10(v) else v
        lo <- if (is.finite(g$min[r])) to_s(g$min[r]) else M[m, "min"]
        hi <- if (is.finite(g$max[r])) to_s(g$max[r]) else M[m, "max"]
        bounds[[m]] <- c(max(M[m, "min"], min(lo, M[m, "max"])),
                         min(M[m, "max"], max(hi, M[m, "min"])))
        on[m] <- isTRUE(as.logical(g$enabled[r]))
      }
    }
    f <- previous$flags
    if (is.data.frame(f) && all(c("flag", "enabled") %in% names(f))) {
      keep <- f$flag %in% names(flags_on)
      flags_on[f$flag[keep]] <- as.logical(f$enabled[keep])
    }
    o <- previous$overrides
    if (length(o) && !is.null(names(o))) {
      o <- o[names(o) %in% qd$ids & o %in% c("keep", "exclude")]
      overrides <- stats::setNames(as.character(o), names(o))
    }
    if (is.numeric(previous$k) && length(previous$k) == 1) k <- previous$k
  }
  list(bounds = bounds, on = on, flags_on = flags_on, overrides = overrides, k = k)
}

# ── Classification ───────────────────────────────────────────────────────────
# A gate pinned to the end of its slider is open on that side, so cells
# beyond the rounded slider range are never failed by accident.
.qc_gate_fail <- function(s, sp, b)
{
  eps <- sp$step / 1000
  lo <- if (b[1] <= sp$min + eps) -Inf else b[1]
  hi <- if (b[2] >= sp$max - eps)  Inf else b[2]
  !is.na(s) & (s < lo - eps | s > hi + eps)
}

.qc_classify <- function(qd, bounds, on, flags_on, overrides)
{
  M <- qd$metrics
  fails <- vapply(M$metric, function(m) {
    if (!isTRUE(on[[m]])) return(logical(qd$n))
    .qc_gate_fail(qd$S[, m], M[m, ], bounds[[m]])
  }, logical(qd$n))
  fails <- matrix(fails, nrow = qd$n, dimnames = list(qd$ids, M$metric))
  ffl <- vapply(qd$flags$flag, function(f) {
    if (!isTRUE(flags_on[[f]])) return(logical(qd$n))
    qd$df[[f]] %in% qd$flags[f, "exclude_if"]
  }, logical(qd$n))
  ffl <- matrix(ffl, nrow = qd$n, dimnames = list(qd$ids, qd$flags$flag))

  gate_pass <- rowSums(fails) == 0 & rowSums(ffl) == 0
  ov <- overrides[qd$ids]
  status <- ifelse(!is.na(ov) & ov == "exclude", "excluded",
            ifelse(!is.na(ov) & ov == "keep" & !gate_pass, "kept",
            ifelse(gate_pass, "pass", "fail")))
  names(status) <- qd$ids
  list(fails = fails, flag_fails = ffl, gate_pass = gate_pass,
       pass = stats::setNames(status %in% c("pass", "kept"), qd$ids), status = status)
}

# Labels of the gates/flags each cell fails, e.g. "MAPD, Spikiness".
.qc_failed_labels <- function(cls, qd)
{
  lab <- c(qd$metrics$label, qd$flags$label)
  all_f <- cbind(cls$fails, cls$flag_fails)
  if (ncol(all_f) == 0) return(rep("", qd$n))
  apply(all_f, 1, function(r) paste(lab[r], collapse = ", "))
}

# Gates as a data frame on the original scale (for export and return value).
.qc_gates_table <- function(qd, bounds, on)
{
  M <- qd$metrics
  to_x <- function(m, v) if (M[m, "log"]) 10^v else v
  data.frame(
    metric = M$metric, label = M$label, enabled = as.logical(on[M$metric]),
    direction = M$dir,
    min = vapply(M$metric, function(m) {
      b <- bounds[[m]][1]; if (b <= M[m, "min"]) -Inf else to_x(m, b)
    }, numeric(1)),
    max = vapply(M$metric, function(m) {
      b <- bounds[[m]][2]; if (b >= M[m, "max"]) Inf else to_x(m, b)
    }, numeric(1)),
    data_min = vapply(M$metric, function(m) min(qd$df[[m]], na.rm = TRUE), numeric(1)),
    data_max = vapply(M$metric, function(m) max(qd$df[[m]], na.rm = TRUE), numeric(1)),
    n_fail = NA_integer_, stringsAsFactors = FALSE, row.names = NULL)
}

.qc_result <- function(qd, state)
{
  cls <- .qc_classify(qd, state$bounds, state$on, state$flags_on, state$overrides)
  gates <- .qc_gates_table(qd, state$bounds, state$on)
  gates$n_fail <- as.integer(colSums(cls$fails))
  flags <- data.frame(flag = qd$flags$flag, label = qd$flags$label,
                      enabled = as.logical(state$flags_on[qd$flags$flag]),
                      exclude_if = qd$flags$exclude_if,
                      n_fail = as.integer(colSums(cls$flag_fails)),
                      stringsAsFactors = FALSE, row.names = NULL)
  ov <- state$overrides[qd$ids]
  cells <- cbind(qd$df["cell_id"],
                 data.frame(qc_pass = unname(cls$pass), status = unname(cls$status),
                            override = ifelse(is.na(ov), "", unname(ov)),
                            failed_metrics = .qc_failed_labels(cls, qd),
                            stringsAsFactors = FALSE),
                 qd$df[setdiff(names(qd$df), "cell_id")])
  rownames(cells) <- qd$ids
  list(pass = qd$ids[cls$pass], fail = qd$ids[!cls$pass], cells = cells,
       thresholds = gates, flags = flags, overrides = state$overrides, k = state$k)
}

# Store a result in res, the way getQCs() stores QC_metrics: res$pass and
# res$fail (cell ids), res$cells (per-cell table) and res$QC_session
# (thresholds, flags and manual decisions, which launch_ui_qc() resumes from).
.qc_attach <- function(res, out)
{
  res$pass  <- out$pass
  res$fail  <- out$fail
  res$cells <- out$cells
  res$QC_session <- out[c("thresholds", "flags", "overrides", "k")]
  res
}

# ── Sparkline (inline SVG, coloured by CSS so theme switches are free) ───────
.qc_spark_svg <- function(h, lo, hi, on = TRUE, marker = NULL, height = 28,
                          class = "qc-spark")
{
  W <- 240
  n  <- length(h$counts)
  bw <- W / n
  mx <- max(h$counts, 1)
  ht <- ifelse(h$counts > 0, pmax(1.5, (height - 2) * h$counts / mx), 0)
  inside <- h$mids >= lo & h$mids <= hi
  cls <- if (!on) rep("off", n) else ifelse(inside, "in", "out")
  bars <- sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="1" class="%s"/>',
                  (seq_len(n) - 1) * bw + 0.6, height - ht, max(bw - 1.2, 0.8), ht, cls)
  mk <- ""
  if (!is.null(marker) && is.finite(marker)) {
    rng <- range(h$breaks)
    mx_x <- W * (min(max(marker, rng[1]), rng[2]) - rng[1]) / diff(rng)
    mk <- sprintf('<line x1="%.2f" x2="%.2f" y1="0" y2="%d" class="mk"/>', mx_x, mx_x, height)
  }
  shiny::HTML(sprintf(
    '<svg class="%s" viewBox="0 0 %d %d" preserveAspectRatio="none" aria-hidden="true">%s<line x1="0" x2="%d" y1="%d" y2="%d" class="base"/>%s</svg>',
    class, W, height, paste(bars, collapse = ""), W, height, height, mk))
}

# Readout of a gate in words: "≤ 0.243", "≥ 65.3k", "0.1 – 0.4", "no limit".
.qc_gate_text <- function(sp, b)
{
  to_x <- function(v) if (sp$log) 10^v else v
  eps <- sp$step / 1000
  lo_open <- b[1] <= sp$min + eps
  hi_open <- b[2] >= sp$max - eps
  if (lo_open && hi_open) return("no limit")
  if (lo_open) return(paste(.qc_s[["le"]], .qc_fmt(to_x(b[2]))))
  if (hi_open) return(paste(.qc_s[["ge"]], .qc_fmt(to_x(b[1]))))
  paste(.qc_fmt(to_x(b[1])), .qc_s[["ndash"]], .qc_fmt(to_x(b[2])))
}

# ── Theme ────────────────────────────────────────────────────────────────────
.qc_font <- "-apple-system, BlinkMacSystemFont, 'SF Pro Text', 'Inter', 'Segoe UI', system-ui, sans-serif"
.qc_mono <- "ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace"

.qc_bs_theme <- function()
{
  bslib::bs_theme(
    version = 5, preset = "bootstrap",
    primary = "#0A84FF", success = "#34C759", danger = "#FF3B30",
    warning = "#FF9500", info = "#0A84FF",
    base_font = bslib::font_collection("-apple-system", "BlinkMacSystemFont",
                                       "SF Pro Text", "Inter", "Segoe UI",
                                       "system-ui", "sans-serif"),
    code_font = bslib::font_collection("ui-monospace", "SF Mono", "SFMono-Regular",
                                       "Menlo", "Consolas", "monospace"),
    "font-size-base" = "0.875rem",
    "border-radius" = "0.6rem",
    "border-radius-lg" = "1rem",
    "border-radius-sm" = "0.4rem"
  )
}

# Copy-number colours from ASCATsc_plot.R. CN 2 and 3 are darkened on light
# backgrounds, as in the ASCAT.sc Explorer web app.
.qc_cn_colors <- function(mode)
{
  cols <- get_total_cn_colors()
  if (mode == "light") { cols["2"] <- "#888888"; cols["3"] <- "#c49631" }
  cols
}

# Plot colours per theme. Green = pass and red = fail are the only strong
# colours; gate shading, the selected cell and magnitude scales are neutral.
# Pass/fail pairs clear the colour-blind separation target on each card
# surface (dataviz validator: dE 9.4 light, 11.2 dark; target >= 8).
.qc_plot_palette <- function(mode = c("dark", "light"))
{
  mode <- match.arg(mode)
  common <- list(mode = mode, font = .qc_font)
  if (mode == "dark")
    c(common, list(
      pass = "#12A37A", fail = "#FF453A",
      surface = "#242426", text = "#f5f5f7", text2 = "rgba(235,235,245,0.78)",
      muted = "rgba(235,235,245,0.55)", grid = "#323235", axis = "#48484a",
      band = "rgba(255,255,255,0.035)", raw = "rgba(235,235,245,0.30)",
      seg = "rgba(235,235,245,0.62)", ring = "#f5f5f7",
      gate_wash = "rgba(235,235,245,0.05)",
      hover_bg = "#2c2c2e", hover_border = "rgba(255,255,255,0.12)",
      seq = c("#636366", "#f2f2f7"), mag = c("#2c2c2e", "#f2f2f7"),
      cn = .qc_cn_colors("dark")))
  else
    c(common, list(
      pass = "#058A5E", fail = "#D92D20",
      surface = "#ffffff", text = "#1d1d1f", text2 = "rgba(60,60,67,0.78)",
      muted = "rgba(60,60,67,0.6)", grid = "#ececf0", axis = "#d1d1d6",
      band = "rgba(0,0,0,0.028)", raw = "rgba(60,60,67,0.28)",
      seg = "rgba(29,29,31,0.55)", ring = "#1d1d1f",
      gate_wash = "rgba(60,60,67,0.045)",
      hover_bg = "#ffffff", hover_border = "rgba(0,0,0,0.10)",
      seq = c("#c7c7cc", "#1d1d1f"), mag = c("#f2f2f7", "#1d1d1f"),
      cn = .qc_cn_colors("light")))
}

.qc_axis <- function(pal, title = NULL, ...)
  utils::modifyList(list(
    title = if (is.null(title)) NULL else list(text = title, font = list(size = 11, color = pal$muted)),
    color = pal$muted, gridcolor = pal$grid, linecolor = pal$axis, zeroline = FALSE,
    showline = TRUE, ticks = "", tickfont = list(size = 10, color = pal$muted),
    automargin = TRUE), list(...))

.qc_layout <- function(p, pal, ...)
  plotly::layout(p,
    font = list(family = pal$font, color = pal$text, size = 11),
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    hoverlabel = list(bgcolor = pal$hover_bg, bordercolor = pal$hover_border,
                      font = list(family = pal$font, color = pal$text, size = 12)),
    ...)

# Log-scale tick labels for metrics gated on log10.
.qc_log_ticks <- function(rng)
{
  e <- seq(floor(rng[1] * 2) / 2, ceiling(rng[2] * 2) / 2, by = 0.5)
  if (length(e) > 7) e <- e[e == round(e)]
  list(tickvals = e, ticktext = .qc_fmt(10^e, 2))
}

# ── Per-cell copy-number profile ─────────────────────────────────────────────
# Bins are put on the copy-number scale exactly as plotSolution() does:
# centred by normaliseByPloidy()'s mean segment value, then converted with
# transform_bulk2tumour(), so they line up with the fitted segments.
.qc_profile_data <- function(res, j, sex = NULL)
{
  if (is.na(j)) return(NULL)
  tr   <- res$allTracks.processed[[j]]
  prof <- res$allProfiles[[j]]
  sol  <- res$allSolutions[[j]]
  if (!is.list(tr) || is.null(tr$lCTS) || !is.data.frame(prof) || !is.list(sol))
    return(NULL)
  rho <- as.numeric(.qc_or(sol$purity, 1))
  psi <- as.numeric(.qc_or(sol$ploidy, 2))
  lcts <- tr$lCTS
  chrs <- .qc_or(names(lcts), as.character(res$chr))
  if (length(chrs) != length(lcts)) return(NULL)

  bins <- do.call(rbind, lapply(seq_along(lcts), function(k) {
    x <- lcts[[k]]
    data.frame(chr = chrs[k], start = as.numeric(x$start), end = as.numeric(x$end),
               smoothed = as.numeric(x$smoothed), stringsAsFactors = FALSE)
  }))
  meanS <- if (!is.null(tr$lSegs)) mean(unlist(lapply(tr$lSegs, function(s)
    rep(s$output$seg.mean, s$output$num.mark))), na.rm = TRUE) else 0
  if (!is.finite(meanS)) meanS <- 0
  ismale <- identical(sex, "male")
  isX <- gsub("chr", "", bins$chr) %in% c("X", "23")
  bins$cn <- transform_bulk2tumour(bins$smoothed - meanS, rho, psi,
                                   ismale = ismale, isX = isX,
                                   isPON = isTRUE(res$isPON))

  segs <- data.frame(chr = as.character(prof$chromosome),
                     start = as.numeric(prof$start), end = as.numeric(prof$end),
                     cn = as.numeric(prof$total_copy_number),
                     cn_raw = as.numeric(prof$total_copy_number_logr),
                     stringsAsFactors = FALSE)

  # Genome coordinates (Mb), chromosomes in track order
  chr_len <- vapply(chrs, function(cc) max(c(bins$end[bins$chr == cc],
                                             segs$end[segs$chr == cc]), na.rm = TRUE), numeric(1))
  offs <- stats::setNames(c(0, cumsum(chr_len)[-length(chr_len)]), chrs)
  bins$x <- (offs[bins$chr] + (bins$start + bins$end) / 2) / 1e6
  segs <- segs[segs$chr %in% chrs, , drop = FALSE]
  segs$x0 <- (offs[segs$chr] + segs$start) / 1e6
  segs$x1 <- (offs[segs$chr] + segs$end) / 1e6

  # Segment of each bin (for one hover layer carrying both values)
  bins$seg_cn <- NA_real_
  for (cc in chrs) {
    bi <- which(bins$chr == cc); si <- which(segs$chr == cc)
    if (!length(bi) || !length(si)) next
    k <- findInterval((bins$start[bi] + bins$end[bi]) / 2, segs$start[si])
    k[k == 0] <- 1
    bins$seg_cn[bi] <- segs$cn[si][k]
  }

  chr_tab <- data.frame(chr = gsub("chr", "", chrs), x0 = offs / 1e6,
                        x1 = (offs + chr_len) / 1e6, stringsAsFactors = FALSE)
  chr_tab$mid <- (chr_tab$x0 + chr_tab$x1) / 2
  ymax <- min(10, max(5, ceiling(max(segs$cn, na.rm = TRUE)) + 1))
  list(bins = bins, segs = segs, chr = chr_tab, purity = rho, ploidy = psi,
       ymax = ymax, sex = sex)
}

.qc_plot_profile <- function(pd, pal, compact = FALSE)
{
  b <- pd$bins
  ymax <- pd$ymax
  yclip <- pmin(pmax(b$cn, 0), ymax + 0.3)
  hover <- sprintf("<b>CN %s</b>  segment<br>%s  bin<br>chr%s:%s Mb",
                   ifelse(is.na(b$seg_cn), .qc_s[["dash"]], b$seg_cn),
                   ifelse(is.finite(b$cn), formatC(b$cn, format = "f", digits = 2), .qc_s[["dash"]]),
                   gsub("chr", "", b$chr),
                   formatC(((b$start + b$end) / 2) / 1e6, format = "f", digits = 1))

  shapes <- lapply(seq_len(nrow(pd$chr)), function(k) {
    if (k %% 2 == 1) return(NULL)
    list(type = "rect", xref = "x", yref = "paper", x0 = pd$chr$x0[k], x1 = pd$chr$x1[k],
         y0 = 0, y1 = 1, fillcolor = pal$band, line = list(width = 0), layer = "below")
  })
  shapes <- Filter(Negate(is.null), shapes)

  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, type = "scattergl", mode = "markers", x = b$x, y = yclip,
                         marker = list(size = if (compact) 2.5 else 3, color = pal$raw),
                         text = hover, hovertemplate = "%{text}<extra></extra>",
                         name = "Bins", showlegend = FALSE)
  s <- pd$segs
  # Non-integer segment means: the gap to the integer bar is the non-integerness
  ok <- is.finite(s$cn_raw)
  if (any(ok)) {
    yr <- pmin(pmax(s$cn_raw[ok], 0), ymax + 0.3)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines",
                           x = c(rbind(s$x0[ok], s$x1[ok], NA)), y = c(rbind(yr, yr, NA)),
                           line = list(color = pal$seg, width = 1.5), hoverinfo = "skip",
                           name = "Segment mean", showlegend = FALSE)
  }
  for (st in sort(unique(s$cn[is.finite(s$cn)]))) {
    k <- which(s$cn == st)
    key <- as.character(min(6, max(0, st)))                  # 6+ share a colour
    col <- .qc_or(unname(pal$cn[key]), pal$seg)
    yv <- min(st, ymax + 0.3)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines",
                           x = c(rbind(s$x0[k], s$x1[k], NA)), y = rep(c(yv, yv, NA), length(k)),
                           line = list(color = col, width = if (compact) 4 else 5),
                           hoverinfo = "skip", name = paste("CN", st), showlegend = FALSE)
  }
  p <- .qc_layout(p, pal,
    showlegend = FALSE, hovermode = "closest", shapes = shapes,
    margin = list(l = 44, r = 10, t = 6, b = 30),
    xaxis = .qc_axis(pal, tickvals = pd$chr$mid, ticktext = pd$chr$chr, showgrid = FALSE,
                     range = c(0, max(pd$chr$x1)), tickfont = list(size = if (compact) 9 else 10,
                     color = pal$muted), fixedrange = FALSE),
    yaxis = .qc_axis(pal, "Copy number", range = c(-0.25, ymax + 0.45), dtick = 1,
                     tick0 = 0, showline = FALSE))
  plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                 modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d",
                                            "hoverClosestCartesian", "hoverCompareCartesian",
                                            "toggleSpikelines"),
                 toImageButtonOptions = list(format = "png", scale = 2))
}

# ── Explore scatter (2D / 3D) ────────────────────────────────────────────────
# Metrics are plotted on the slider scale (log10 for read counts), so gate
# boxes and axes share one coordinate system; hovers show real values.
.qc_plot_scatter <- function(qd, xm, ym, zm = NULL, status, pal, bounds, on,
                             selected = NULL, color_by = "status", dim = "2d")
{
  M  <- qd$metrics
  ms <- c(xm, ym, if (dim == "3d") zm)
  raw <- function(m) qd$df[[m]]
  hover <- sprintf("<b>%s</b><br>%s", qd$ids, do.call(paste, c(lapply(ms, function(m)
    sprintf("%s  %s", .qc_fmt(raw(m)), M[m, "label"])), sep = "<br>")))
  n  <- qd$n
  sz <- if (n <= 1000) 8 else if (n <= 5000) 6 else 4.5
  # 2D opens on the slider window when extreme cells lie beyond it (they are
  # counted in an annotation; Autoscale shows them).
  clipped <- function(m) M[m, "n_above"] + M[m, "n_below"] > 0
  ax <- function(m) {
    sp <- M[m, ]
    a <- .qc_axis(pal, sp$label)
    win <- if (dim == "2d" && clipped(m)) c(sp$min, sp$max) else range(qd$S[, m], na.rm = TRUE)
    if (dim == "2d" && clipped(m)) a$range <- win + c(-1, 1) * 0.04 * diff(win)
    if (sp$log) a <- c(a, .qc_log_ticks(win))
    a
  }
  data_rng <- function(m) range(qd$S[, m], na.rm = TRUE)
  # Gate window on the slider scale; open sides extend to the data range
  gate_win <- function(m) {
    sp <- M[m, ]; r <- data_rng(m)
    if (!isTRUE(on[[m]])) return(r)
    b <- bounds[[m]]; eps <- sp$step / 1000
    c(if (b[1] <= sp$min + eps) min(r[1], b[1]) else b[1],
      if (b[2] >= sp$max - eps) max(r[2], b[2]) else b[2])
  }
  groups <- list(
    list(key = "pass",     name = "Pass",              col = pal$pass, sym = "circle"),
    list(key = "fail",     name = "Fail",              col = pal$fail, sym = "circle"),
    list(key = "kept",     name = "Kept (manual)",     col = pal$pass, sym = "diamond"),
    list(key = "excluded", name = "Excluded (manual)", col = pal$fail, sym = "diamond"))

  p <- plotly::plot_ly(source = "explore")
  type <- if (dim == "3d") "scatter3d" else "scattergl"
  add_pts <- function(p, k, ...) {
    args <- list(p, type = type, mode = "markers",
                 x = qd$S[k, xm], y = qd$S[k, ym], customdata = qd$ids[k],
                 text = hover[k], hovertemplate = "%{text}<extra></extra>", ...)
    if (dim == "3d") args$z <- qd$S[k, zm]
    do.call(plotly::add_trace, args)
  }
  if (color_by == "status") {
    for (g in groups) {
      k <- which(status == g$key)
      if (!length(k)) next
      p <- add_pts(p, k, name = sprintf("%s  %s", g$name, format(length(k), big.mark = ",")),
                   marker = list(size = if (dim == "3d") sz * 0.55 else sz, color = g$col,
                                 symbol = g$sym, opacity = if (g$key %in% c("fail", "excluded")) 0.95 else 0.7,
                                 line = list(color = pal$surface, width = if (dim == "3d") 0 else 1)))
    }
  } else {
    cv <- qd$S[, color_by]
    p <- add_pts(p, seq_len(n), name = M[color_by, "label"], showlegend = FALSE,
                 marker = list(size = if (dim == "3d") sz * 0.55 else sz, color = cv,
                               colorscale = list(c(0, pal$seq[1]), c(1, pal$seq[2])),
                               opacity = 0.85,
                               line = list(color = pal$surface, width = if (dim == "3d") 0 else 1),
                               colorbar = c(list(title = list(text = M[color_by, "label"],
                                                              font = list(size = 11, color = pal$muted)),
                                                 thickness = 10, len = 0.6, outlinewidth = 0,
                                                 tickfont = list(size = 10, color = pal$muted)),
                                            if (M[color_by, "log"])
                                              .qc_log_ticks(range(cv, na.rm = TRUE)))))
  }
  if (!is.null(selected) && selected %in% qd$ids) {
    k <- match(selected, qd$ids)
    # a neutral ring around the cell, so its own pass/fail colour still shows
    p <- add_pts(p, k, name = "Selected", showlegend = FALSE,
                 marker = if (dim == "3d")
                   list(size = 11, symbol = "circle-open", color = pal$ring, opacity = 1)
                 else
                   list(size = 17, color = "rgba(0,0,0,0)", opacity = 1,
                        line = list(color = pal$ring, width = 2.5)))
  }

  # 3D: the gates on the three axes as a see-through (neutral) box
  if (dim == "3d" && any(vapply(c(xm, ym, zm), function(m) isTRUE(on[[m]]), logical(1)))) {
    gx <- gate_win(xm); gy <- gate_win(ym); gz <- gate_win(zm)
    vx <- gx[c(1, 2, 2, 1, 1, 2, 2, 1)]
    vy <- gy[c(1, 1, 2, 2, 1, 1, 2, 2)]
    vz <- gz[c(1, 1, 1, 1, 2, 2, 2, 2)]
    p <- plotly::add_trace(p, type = "mesh3d", x = vx, y = vy, z = vz,
                           i = c(0, 0, 4, 4, 0, 0, 3, 3, 0, 0, 1, 1),
                           j = c(1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 2, 6),
                           k = c(2, 3, 6, 7, 5, 4, 6, 7, 7, 4, 6, 5),
                           # facecolor, not color: plotly's R wrapper remaps `color`
                           facecolor = rep(pal$ring, 12), opacity = 0.05, flatshading = TRUE,
                           hoverinfo = "skip", showlegend = FALSE, showscale = FALSE,
                           name = "Thresholds")
    e <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(5, 6), c(6, 7), c(7, 8), c(8, 5),
               c(1, 5), c(2, 6), c(3, 7), c(4, 8))
    seg <- function(v) c(t(cbind(v[e[, 1]], v[e[, 2]], NA)))
    p <- plotly::add_trace(p, type = "scatter3d", mode = "lines", x = seg(vx), y = seg(vy),
                           z = seg(vz), line = list(color = pal$muted, width = 2, dash = "dash"),
                           hoverinfo = "skip", showlegend = TRUE, name = "Thresholds")
  }

  uirev <- paste(dim, xm, ym, zm, color_by)
  leg <- list(orientation = "h", x = 0, xanchor = "left", y = 1.02, yanchor = "bottom",
              itemsizing = "constant", font = list(size = 12, color = pal$text2),
              bgcolor = "rgba(0,0,0,0)")

  if (dim == "3d") {
    sax <- function(m) {
      a <- list(title = list(text = M[m, "label"], font = list(size = 11, color = pal$muted)),
                color = pal$muted, gridcolor = pal$grid, zerolinecolor = pal$grid,
                backgroundcolor = "rgba(0,0,0,0)", showbackground = FALSE,
                tickfont = list(size = 10, color = pal$muted))
      if (M[m, "log"]) a <- c(a, .qc_log_ticks(range(qd$S[, m], na.rm = TRUE)))
      a
    }
    p <- .qc_layout(p, pal, legend = leg, uirevision = uirev,
                    margin = list(l = 0, r = 0, t = 30, b = 0),
                    scene = list(xaxis = sax(xm), yaxis = sax(ym), zaxis = sax(zm),
                                 camera = list(eye = list(x = 1.3, y = 1.3, z = 0.85))))
  } else {
    shapes <- list()
    box <- c(gate_win(xm), gate_win(ym))
    if (isTRUE(on[[xm]]) || isTRUE(on[[ym]]))
      shapes[[1]] <- list(type = "rect", xref = "x", yref = "y", x0 = box[1], x1 = box[2],
                          y0 = box[3], y1 = box[4], fillcolor = pal$gate_wash,
                          line = list(width = 0), layer = "below")
    thr <- function(m, axis) {
      if (!isTRUE(on[[m]])) return(list())
      b <- bounds[[m]]; sp <- M[m, ]; out <- list()
      for (v in b[c(b[1] > sp$min + sp$step / 1000, b[2] < sp$max - sp$step / 1000)])
        out[[length(out) + 1]] <- if (axis == "x")
          list(type = "line", xref = "x", yref = "paper", x0 = v, x1 = v, y0 = 0, y1 = 1,
               line = list(color = pal$muted, width = 1, dash = "dash"), layer = "below")
        else
          list(type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1, y0 = v, y1 = v,
               line = list(color = pal$muted, width = 1, dash = "dash"), layer = "below")
      out
    }
    shapes <- c(shapes, thr(xm, "x"), thr(ym, "y"))
    ann <- list()
    note <- function(x, y, xa, ya, txt)
      list(xref = "paper", yref = "paper", x = x, y = y, xanchor = xa, yanchor = ya,
           showarrow = FALSE, text = txt, font = list(size = 11, color = pal$muted))
    fc <- function(v) format(v, big.mark = ",")
    if (M[ym, "n_above"]) ann <- c(ann, list(note(0.01, 0.99, "left", "top",
                                                   sprintf(paste0(.qc_s[["uarr"]], " %s more cells above"), fc(M[ym, "n_above"])))))
    if (M[ym, "n_below"]) ann <- c(ann, list(note(0.01, 0.01, "left", "bottom",
                                                   sprintf(paste0(.qc_s[["darr"]], " %s more cells below"), fc(M[ym, "n_below"])))))
    if (M[xm, "n_above"]) ann <- c(ann, list(note(0.99, 0.01, "right", "bottom",
                                                   sprintf(paste0("%s more cells to the right ", .qc_s[["rarr"]]), fc(M[xm, "n_above"])))))
    if (M[xm, "n_below"]) ann <- c(ann, list(note(0.01, 0.06, "left", "bottom",
                                                   sprintf(paste0(.qc_s[["larr"]], " %s more cells to the left"), fc(M[xm, "n_below"])))))
    if (!is.null(selected) && selected %in% qd$ids) {
      k <- match(selected, qd$ids)
      if (is.finite(qd$S[k, xm]) && is.finite(qd$S[k, ym]))
        ann[[length(ann) + 1]] <- list(x = qd$S[k, xm], y = qd$S[k, ym], text = selected, showarrow = TRUE,
                         arrowhead = 0, arrowwidth = 1, arrowcolor = pal$muted, ax = 18, ay = -22,
                         font = list(size = 11, color = pal$text), bgcolor = pal$hover_bg,
                         bordercolor = pal$hover_border, borderpad = 3)
    }
    p <- .qc_layout(p, pal, legend = leg, uirevision = uirev, shapes = shapes,
                    annotations = ann, hovermode = "closest", hoverdistance = 24,
                    dragmode = "zoom", margin = list(l = 10, r = 10, t = 36, b = 10),
                    xaxis = ax(xm), yaxis = ax(ym))
  }
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                      modeBarButtonsToRemove = c("hoverClosestCartesian", "hoverCompareCartesian",
                                                 "toggleSpikelines"),
                      toImageButtonOptions = list(format = "png", scale = 2))
  p <- plotly::event_register(p, "plotly_click")
  p <- plotly::event_register(p, "plotly_selected")
  plotly::event_register(p, "plotly_deselect")
}

# ── Distribution of one metric, stacked by overall pass / fail ───────────────
.qc_plot_hist <- function(qd, m, b, on, pass, pal)
{
  sp <- qd$metrics[m, ]
  s  <- qd$S[, m]
  br <- seq(sp$min, sp$max, length.out = 41)
  cp <- .qc_hist_counts(s[pass],  br)
  cf <- .qc_hist_counts(s[!pass], br)
  mids <- (br[-1] + br[-length(br)]) / 2
  w  <- br[2] - br[1]
  to_x <- function(v) if (sp$log) 10^v else v
  rng_txt <- sprintf(paste0("%s ", .qc_s[["ndash"]], " %s"), .qc_fmt(to_x(br[-length(br)])), .qc_fmt(to_x(br[-1])))

  p <- plotly::plot_ly()
  p <- plotly::add_bars(p, x = mids, y = cp, width = w, name = "Pass",
                        marker = list(color = pal$pass, line = list(color = pal$surface, width = 1)),
                        text = rng_txt, textposition = "none", hovertemplate = "<b>%{y}</b> pass<br>%{text}<extra></extra>")
  p <- plotly::add_bars(p, x = mids, y = cf, width = w, name = "Fail",
                        marker = list(color = pal$fail, line = list(color = pal$surface, width = 1)),
                        text = rng_txt, textposition = "none", hovertemplate = "<b>%{y}</b> fail<br>%{text}<extra></extra>")
  shapes <- list()
  eps <- sp$step / 1000
  if (isTRUE(on)) {
    if (b[1] > sp$min + eps)
      shapes <- c(shapes, list(
        list(type = "rect", xref = "x", yref = "paper", x0 = sp$min, x1 = b[1], y0 = 0, y1 = 1,
             fillcolor = pal$gate_wash, line = list(width = 0), layer = "below"),
        list(type = "line", xref = "x", yref = "paper", x0 = b[1], x1 = b[1], y0 = 0, y1 = 1,
             line = list(color = pal$text2, width = 1, dash = "dash"))))
    if (b[2] < sp$max - eps)
      shapes <- c(shapes, list(
        list(type = "rect", xref = "x", yref = "paper", x0 = b[2], x1 = sp$max, y0 = 0, y1 = 1,
             fillcolor = pal$gate_wash, line = list(width = 0), layer = "below"),
        list(type = "line", xref = "x", yref = "paper", x0 = b[2], x1 = b[2], y0 = 0, y1 = 1,
             line = list(color = pal$text2, width = 1, dash = "dash"))))
  }
  xa <- .qc_axis(pal, NULL, range = c(sp$min, sp$max), showgrid = FALSE, fixedrange = TRUE,
                 nticks = 5)
  if (sp$log) xa <- c(xa, .qc_log_ticks(c(sp$min, sp$max)))
  # Cells beyond the window are drawn in the end bins; say how many
  ann <- list()
  edge <- function(x, xa_, txt)
    list(xref = "paper", yref = "paper", x = x, y = 1, xanchor = xa_, yanchor = "top",
         showarrow = FALSE, text = txt, font = list(size = 10, color = pal$muted),
         bgcolor = pal$surface)
  if (sp$n_above) ann <- c(ann, list(edge(1, "right", sprintf(paste0("%d beyond %s ", .qc_s[["rsaquo"]]), sp$n_above,
                                                              .qc_fmt(to_x(sp$max))))))
  if (sp$n_below) ann <- c(ann, list(edge(0, "left", sprintf(paste0(.qc_s[["lsaquo"]], " %d below %s"), sp$n_below,
                                                             .qc_fmt(to_x(sp$min))))))
  p <- .qc_layout(p, pal, barmode = "stack", bargap = 0, showlegend = FALSE, shapes = shapes,
                  annotations = ann,
                  margin = list(l = 30, r = 6, t = 4, b = 22), hovermode = "closest",
                  xaxis = xa, yaxis = .qc_axis(pal, NULL, showline = FALSE, fixedrange = TRUE,
                                               nticks = 4))
  plotly::config(p, displayModeBar = FALSE, responsive = TRUE)
}

# ── Spearman correlation between gate metrics ────────────────────────────────
.qc_plot_corr <- function(qd, pal)
{
  M <- qd$metrics
  R <- suppressWarnings(stats::cor(qd$S, method = "spearman", use = "pairwise.complete.obs"))
  lab <- M$label
  p <- plotly::plot_ly(source = "corr")
  # Neutral strength scale (|rho|): red and green stay reserved for fail / pass.
  # The sign is in the tooltip.
  rho_txt <- matrix(formatC(R, format = "f", digits = 2), nrow(R))
  p <- plotly::add_heatmap(p, z = abs(R), x = lab, y = lab, text = rho_txt,
                           zmin = 0, zmax = 1, xgap = 2, ygap = 2,
                           colorscale = list(c(0, pal$mag[1]), c(1, pal$mag[2])),
                           colorbar = list(thickness = 8, len = 0.7, outlinewidth = 0,
                                           tickvals = c(0, 0.5, 1),
                                           tickfont = list(size = 10, color = pal$muted),
                                           title = list(text = paste0("|", .qc_s[["rho"]], "|"),
                                                        font = list(color = pal$muted))),
                           hovertemplate = paste0("<b>", .qc_s[["rho"]], " = %{text}</b><br>%{y}<br>%{x}<extra></extra>"))
  p <- .qc_layout(p, pal, margin = list(l = 10, r = 10, t = 6, b = 10),
                  xaxis = .qc_axis(pal, NULL, showgrid = FALSE, showline = FALSE, tickangle = -40,
                                   fixedrange = TRUE),
                  yaxis = .qc_axis(pal, NULL, showgrid = FALSE, showline = FALSE,
                                   autorange = "reversed", fixedrange = TRUE))
  p <- plotly::config(p, displayModeBar = FALSE, responsive = TRUE)
  plotly::event_register(p, "plotly_click")
}

# Copy-number colour key for the profile plots (colours come from the CSS, so
# they follow the theme)
.qc_cn_legend <- function()
  shiny::div(class = "qc-cn-legend", shiny::span(class = "qc-cn-title", "Copy number"),
             lapply(c(0:5, "6+"), function(k)
               shiny::span(shiny::tags$i(style = sprintf("background: var(--qc-cn-%s)", sub("+", "", k, fixed = TRUE))),
                           k)))

# ── ASCAT.sc hex logo ────────────────────────────────────────────────────────
# The image is set in the stylesheet: the dark hex in dark mode, the light hex in light
# mode, so a theme switch swaps it without a re-render.
.qc_logo <- function(class = NULL)
  shiny::span(class = paste(c("qc-logo-hex", class), collapse = " "),
              role = "img", `aria-label` = "ASCAT.sc")

# =============================================================================
# Embedded assets: stylesheet, client script and logos. They live here, not in
# separate files, so that this one file is the whole app.
# =============================================================================

.qc_logo_uri <- function(which = c("dark", "light"))
{
  b64 <- switch(match.arg(which), dark = .qc_logo_dark_b64, light = .qc_logo_light_b64)
  paste0("data:image/webp;base64,", gsub("[[:space:]]", "", b64))
}

.qc_css <- function()
{
  css <- gsub("@LOGO_DARK@", .qc_logo_uri("dark"), .qc_css_src, fixed = TRUE)
  gsub("@LOGO_LIGHT@", .qc_logo_uri("light"), css, fixed = TRUE)
}

# ── Stylesheet ───────────────────────────────────────────────────────────────
.qc_css_src <- r"---(
/* =============================================================================
   ASCAT.sc QC app (launch_ui_qc) - design system
   Same tokens as the ASCAT.sc Explorer web app: Apple system colours,
   glass chrome, flat opaque content cards. Every colour is a variable, so
   switching theme (data-bs-theme on <html>) never needs a re-render.
   ============================================================================= */

/* -- Tokens --------------------------------------------------------------- */
html[data-bs-theme="light"], html:not([data-bs-theme]) {
  color-scheme: light;
  --qc-bg: #f5f5f7;
  --qc-card: #ffffff;
  --qc-inset: #f5f5f7;
  --qc-item: #e8e8ed;
  --qc-seg-on: #ffffff;
  --qc-border: rgba(0, 0, 0, 0.08);
  --qc-border-2: rgba(0, 0, 0, 0.14);
  --qc-text: #1d1d1f;
  --qc-text-2: rgba(60, 60, 67, 0.78);
  --qc-text-3: rgba(60, 60, 67, 0.6);
  --qc-text-4: rgba(60, 60, 67, 0.4);
  --qc-track: #dcdce1;
  --qc-switch-off: #d1d1d6;
  /* Colour carries meaning only: green = pass, red = fail. Controls are
     neutral graphite. Pass/fail were checked for colour-blind separation. */
  --qc-pass: #058a5e;
  --qc-fail: #d92d20;
  --qc-good-text: #046c4a;
  --qc-bad-text: #b42318;
  --qc-pass-soft: rgba(5, 138, 94, 0.12);
  --qc-fail-soft: rgba(217, 45, 32, 0.10);
  --qc-fail-mid: rgba(217, 45, 32, 0.32);
  --qc-pass-line: rgba(5, 138, 94, 0.5);
  --qc-fail-line: rgba(217, 45, 32, 0.45);
  --qc-spark-in: rgba(60, 60, 67, 0.42);
  --qc-accent: #48484a;
  --qc-accent-soft: rgba(60, 60, 67, 0.08);
  --qc-accent-text: #1d1d1f;
  --qc-accent-rgb: 60, 60, 67;
  --qc-switch-on: #1d1d1f;
  --qc-primary-bg: #1d1d1f;
  --qc-primary-text: #ffffff;
  --qc-hover: rgba(0, 0, 0, 0.035);
  --qc-focus: rgba(60, 60, 67, 0.45);
  --qc-cn-2: #888888;
  --qc-cn-3: #c49631;
  --qc-glass: rgba(255, 255, 255, 0.72);
  --qc-glass-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.8), 0 4px 24px rgba(0, 0, 0, 0.05);
  --qc-menu: rgba(255, 255, 255, 0.94);
  --qc-menu-shadow: 0 10px 30px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.08);
  --qc-scroll: #c7c7cc;
}
html[data-bs-theme="dark"] {
  color-scheme: dark;
  --qc-bg: #1c1c1e;
  --qc-card: #242426;
  --qc-inset: #1c1c1e;
  --qc-item: #3a3a3c;
  --qc-seg-on: #5a5a5e;
  --qc-border: rgba(255, 255, 255, 0.08);
  --qc-border-2: rgba(255, 255, 255, 0.14);
  --qc-text: #f5f5f7;
  --qc-text-2: rgba(235, 235, 245, 0.78);
  --qc-text-3: rgba(235, 235, 245, 0.55);
  --qc-text-4: rgba(235, 235, 245, 0.35);
  --qc-track: #48484a;
  --qc-switch-off: #39393d;
  /* Colour carries meaning only: green = pass, red = fail. Controls are
     space grey. Pass/fail were checked for colour-blind separation. */
  --qc-pass: #12a37a;
  --qc-fail: #ff453a;
  --qc-good-text: #3fd19f;
  --qc-bad-text: #ff7b72;
  --qc-pass-soft: rgba(18, 163, 122, 0.18);
  --qc-fail-soft: rgba(255, 69, 58, 0.16);
  --qc-fail-mid: rgba(255, 69, 58, 0.36);
  --qc-pass-line: rgba(18, 163, 122, 0.6);
  --qc-fail-line: rgba(255, 69, 58, 0.5);
  --qc-spark-in: rgba(235, 235, 245, 0.40);
  --qc-accent: #aeaeb2;
  --qc-accent-soft: rgba(174, 174, 178, 0.16);
  --qc-accent-text: #f5f5f7;
  --qc-accent-rgb: 174, 174, 178;
  --qc-switch-on: #e5e5ea;
  --qc-primary-bg: #e5e5ea;
  --qc-primary-text: #1c1c1e;
  --qc-hover: rgba(255, 255, 255, 0.05);
  --qc-focus: rgba(235, 235, 245, 0.45);
  --qc-cn-2: #c1c1c1;
  --qc-cn-3: #e3b55f;
  --qc-glass: rgba(28, 28, 30, 0.66);
  --qc-glass-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.06), 0 4px 24px rgba(0, 0, 0, 0.35);
  --qc-menu: rgba(44, 44, 46, 0.96);
  --qc-menu-shadow: 0 10px 30px rgba(0, 0, 0, 0.5), 0 0 0 0.5px rgba(255, 255, 255, 0.1);
  --qc-scroll: #48484a;
}
html[data-bs-theme] {
  /* copy-number colours of the profile plot (get_total_cn_colors, capped at 6+) */
  --qc-cn-0: #153570;
  --qc-cn-1: #577aba;
  --qc-cn-4: #d6804f;
  --qc-cn-5: #b3402e;
  --qc-cn-6: #821010;
  --qc-header-h: 52px;
  --qc-sidebar-w: 324px;
  --qc-mono: ui-monospace, "SF Mono", SFMono-Regular, Menlo, Consolas, monospace;
  /* Bootstrap colour-mode variables, so Bootstrap/DT/selectize parts follow */
  --bs-body-bg: var(--qc-bg);
  --bs-body-color: var(--qc-text);
  --bs-emphasis-color: var(--qc-text);
  --bs-secondary-color: var(--qc-text-3);
  --bs-tertiary-color: var(--qc-text-4);
  --bs-secondary-bg: var(--qc-item);
  --bs-tertiary-bg: var(--qc-inset);
  --bs-border-color: var(--qc-border);
  --bs-border-color-translucent: var(--qc-border);
  --bs-link-color: var(--qc-accent);
  --bs-link-hover-color: var(--qc-accent);
  --bs-link-color-rgb: var(--qc-accent-rgb);
  --bs-primary-rgb: var(--qc-accent-rgb);
  --bs-primary: var(--qc-accent);
}

/* -- Base ----------------------------------------------------------------- */
html, body { background: var(--qc-bg); }
body {
  color: var(--qc-text);
  -webkit-font-smoothing: antialiased;
  letter-spacing: -0.003em;
}
.qc-app { min-height: 100vh; }
.qc-app :focus-visible { outline: 2px solid var(--qc-focus); outline-offset: 2px; }
.qc-mono { font-family: var(--qc-mono); font-variant-numeric: tabular-nums; }
.qc-muted { color: var(--qc-text-3); }
.qc-faint { color: var(--qc-text-4); }

::-webkit-scrollbar { width: 9px; height: 9px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb {
  background: var(--qc-scroll); border-radius: 9px;
  border: 2px solid transparent; background-clip: content-box;
}
* { scrollbar-color: var(--qc-scroll) transparent; scrollbar-width: thin; }

/* Quick updates should not flash; slow widgets dim after a short delay. */
.qc-app .recalculating { opacity: 1; }
.qc-app .html-widget-output.recalculating {
  opacity: 0.6; transition: opacity 0.25s ease 0.25s;
}

/* -- Header (glass chrome, sticky over scrolling content) ----------------- */
.qc-header {
  position: sticky; top: 0; z-index: 1030;
  height: var(--qc-header-h);
  display: flex; align-items: center; gap: 12px;
  padding: 0 16px;
  background: var(--qc-glass);
  -webkit-backdrop-filter: blur(24px) saturate(180%);
  backdrop-filter: blur(24px) saturate(180%);
  border-bottom: 1px solid var(--qc-border);
  box-shadow: var(--qc-glass-shadow);
}
.qc-brand { display: flex; align-items: center; gap: 10px; min-width: 0; }
/* Hex logo: dark hex in dark mode, light hex in light mode (embedded WebP) */
.qc-logo-hex {
  display: inline-block; flex: none; width: 35px; height: 40px;
  background: url("@LOGO_LIGHT@") center / contain no-repeat;
}
html[data-bs-theme="dark"] .qc-logo-hex { background-image: url("@LOGO_DARK@"); }
.qc-title { font-size: 14px; font-weight: 600; letter-spacing: -0.01em; white-space: nowrap; }
.qc-title span { color: var(--qc-text-3); font-weight: 400; }
.qc-pills { display: flex; align-items: center; gap: 6px; flex-wrap: nowrap; overflow: hidden; }
.qc-pill {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 11.5px; line-height: 1; padding: 5px 9px; border-radius: 999px;
  background: var(--qc-item); color: var(--qc-text-3); white-space: nowrap;
}
.qc-pill .qc-dot { width: 6px; height: 6px; }
.qc-header-actions { margin-left: auto; display: flex; align-items: center; gap: 8px; }

.qc-btn, .qc-app .btn.qc-btn {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 12px; font-weight: 500; line-height: 1;
  padding: 7px 11px; border-radius: 8px; border: 0;
  background: var(--qc-item); color: var(--qc-text-2);
  transition: background-color 0.15s ease, color 0.15s ease, filter 0.15s ease;
  white-space: nowrap;
}
.qc-btn:hover, .qc-app .btn.qc-btn:hover { color: var(--qc-text); filter: brightness(1.07); background: var(--qc-item); }
.qc-btn:active { filter: brightness(0.95); }
.qc-app .btn.qc-btn-primary, .qc-btn-primary { background: var(--qc-primary-bg); color: var(--qc-primary-text); }
.qc-app .btn.qc-btn-primary:hover, .qc-btn-primary:hover { background: var(--qc-primary-bg); color: var(--qc-primary-text); filter: brightness(1.06); }
.qc-btn-ghost { background: transparent; }
.qc-btn-sm { padding: 5px 9px; font-size: 11.5px; }
.qc-btn-icon { padding: 7px 8px; }
.qc-btn-danger { color: var(--qc-bad-text); }
.qc-app .dropdown-toggle::after { display: none; }

bslib-input-dark-mode {
  --text-1: var(--qc-text-2) !important;
  --text-2: var(--qc-text-3) !important;
  display: inline-flex; align-items: center; justify-content: center;
  width: 31px; height: 27px; border-radius: 8px; background: var(--qc-item);
}

/* Glass dropdown menus (export) */
.qc-app .dropdown-menu {
  --bs-dropdown-bg: var(--qc-menu);
  background: var(--qc-menu);
  border: 1px solid var(--qc-border); border-radius: 12px;
  box-shadow: var(--qc-menu-shadow);
  padding: 5px; min-width: 230px; font-size: 12.5px;
}
.qc-app .dropdown-header {
  font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.07em;
  color: var(--qc-text-4); padding: 6px 10px 4px;
}
.qc-app .dropdown-item {
  display: flex; align-items: center; gap: 9px;
  border-radius: 8px; padding: 7px 10px; color: var(--qc-text-2);
}
.qc-app .dropdown-item .fa, .qc-app .dropdown-item svg { width: 14px; color: var(--qc-text-3); }
.qc-app .dropdown-item:hover, .qc-app .dropdown-item:focus {
  background: var(--qc-accent-soft); color: var(--qc-text);
}
.qc-app .dropdown-item small { margin-left: auto; color: var(--qc-text-4); }

/* -- Layout --------------------------------------------------------------- */
.qc-body {
  display: grid;
  grid-template-columns: var(--qc-sidebar-w) minmax(0, 1fr);
  transition: grid-template-columns 0.22s ease;
}
body.qc-sidebar-hidden .qc-body { grid-template-columns: 0 minmax(0, 1fr); }
.qc-sidebar {
  position: sticky; top: var(--qc-header-h); align-self: start;
  height: calc(100vh - var(--qc-header-h));
  overflow-y: auto; overflow-x: hidden;
  padding: 14px 12px 28px 14px;
  border-right: 1px solid var(--qc-border);
  background: var(--qc-glass);
}
body.qc-sidebar-hidden .qc-sidebar { visibility: hidden; padding: 0; border: 0; }
.qc-main { position: relative; min-width: 0; padding: 16px 22px 48px; }

/* -- Section labels ------------------------------------------------------- */
.qc-label {
  font-size: 10.5px; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.08em; color: var(--qc-text-3);
}
.qc-side-section { margin-bottom: 14px; }
.qc-side-section > .qc-label { display: flex; align-items: center; justify-content: space-between; margin: 0 2px 8px; }
.qc-group-label { margin: 14px 2px 7px; color: var(--qc-text-4); }

/* Auto-threshold row */
.qc-auto {
  display: flex; align-items: center; gap: 6px; flex-wrap: wrap;
  padding: 9px 10px; border-radius: 12px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  font-size: 12px; color: var(--qc-text-2);
}
.qc-auto .form-group { margin: 0; width: auto; }
.qc-auto input.form-control {
  width: 52px; height: 26px; padding: 2px 6px; text-align: center;
  font-family: var(--qc-mono); font-size: 12px;
  background: var(--qc-item); border: 0; border-radius: 7px; color: var(--qc-text);
}
.qc-auto .qc-auto-actions { margin-left: auto; display: flex; gap: 6px; }

/* -- Gate cards ----------------------------------------------------------- */
.qc-gate {
  position: relative;
  padding: 9px 10px 4px; margin-bottom: 7px;
  border-radius: 12px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  transition: border-color 0.15s ease;
}
.qc-gate:hover { border-color: var(--qc-border-2); }
.qc-gate-head {
  display: flex; align-items: center; gap: 8px;
  min-height: 20px; padding-right: 64px;
}
.qc-gate-head .shiny-input-container { width: auto !important; margin: 0; }
.qc-gate-name { font-size: 12.5px; font-weight: 500; color: var(--qc-text); line-height: 1.2; }
.qc-help { color: var(--qc-text-4); font-size: 11px; cursor: help; }
.qc-help:hover { color: var(--qc-text-2); }
.qc-gate-badge {
  position: absolute; top: 9px; right: 10px;
  font-size: 11px; line-height: 1; padding: 4px 7px; border-radius: 999px;
  font-variant-numeric: tabular-nums;
  display: inline-flex; align-items: center; gap: 5px;
  background: var(--qc-item); color: var(--qc-text-2);
}
.qc-gate-badge::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: var(--qc-fail); }
.qc-gate-badge.is-zero { color: var(--qc-text-4); }
.qc-gate-badge.is-zero::before { background: var(--qc-text-4); }
.qc-gate-read {
  display: flex; justify-content: space-between; align-items: baseline;
  margin: 6px 1px 3px; font-size: 11.5px;
}
.qc-gate-limit { font-family: var(--qc-mono); font-variant-numeric: tabular-nums; color: var(--qc-text-2); }
.qc-gate-med { font-size: 10.5px; color: var(--qc-text-4); }
.qc-gate .shiny-input-container:has(.irs) { margin: 0 0 2px; }
.qc-gate:has(.form-check-input:not(:checked)) .qc-gate-body { opacity: 0.42; }
.qc-gate:has(.form-check-input:not(:checked)) .irs--shiny .irs-bar { background: var(--qc-text-4); }
.qc-gate:has(.form-check-input:not(:checked)) .qc-gate-badge { display: none; }

/* Sparklines: neutral inside the gate, red for the cells it removes. */
.qc-spark { display: block; width: 100%; height: 28px; overflow: visible; }
.qc-spark rect.in { fill: var(--qc-spark-in); }
.qc-spark rect.out { fill: var(--qc-fail); opacity: 0.85; }
.qc-spark rect.off { fill: var(--qc-text-4); }
.qc-spark line.base { stroke: var(--qc-border-2); stroke-width: 1; vector-effect: non-scaling-stroke; }
.qc-spark line.mk { stroke: var(--qc-text); stroke-width: 2.5; vector-effect: non-scaling-stroke; }

/* Flags */
.qc-flag {
  display: flex; align-items: center; gap: 9px;
  padding: 8px 10px; margin-bottom: 6px; border-radius: 10px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  font-size: 12.5px; color: var(--qc-text-2);
}
.qc-flag .shiny-input-container { width: auto !important; margin: 0; }
.qc-flag-n { margin-left: auto; font-size: 11px; color: var(--qc-text-4); font-variant-numeric: tabular-nums; }

.qc-side-foot { display: flex; gap: 6px; margin-top: 12px; }
.qc-side-foot .btn { flex: 1; justify-content: center; }
.qc-kbd-hints { margin-top: 16px; font-size: 11px; color: var(--qc-text-4); line-height: 1.9; }
kbd, .qc-kbd {
  font-family: var(--qc-mono); font-size: 10.5px; padding: 1px 5px; border-radius: 5px;
  background: var(--qc-item); color: var(--qc-text-3); box-shadow: none;
}

/* -- iOS switch (bslib input_switch) -------------------------------------- */
.qc-app .form-switch { padding-left: 0; margin: 0; min-height: 0; display: flex; align-items: center; }
.qc-app .form-switch .form-check-input,
.qc-app .form-switch .form-check-input:focus,
.qc-app .form-switch .form-check-input:checked {
  --bs-form-switch-bg: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3e%3ccircle r='3.25' fill='%23fff'/%3e%3c/svg%3e");
}
.qc-app .form-switch .form-check-input {
  width: 30px; height: 18px; margin: 0; flex: none; float: none;
  border: 0; border-radius: 9px; cursor: pointer;
  background-color: var(--qc-switch-off);
  box-shadow: none;
  transition: background-position 0.18s ease, background-color 0.18s ease;
}
.qc-app .form-switch .form-check-input:checked { background-color: var(--qc-switch-on); }
/* Dark mode: the "on" track is bright, so its knob is dark, and on/off can never be confused */
html[data-bs-theme="dark"] .qc-app .form-switch .form-check-input:checked {
  --bs-form-switch-bg: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3e%3ccircle r='3.25' fill='%231c1c1e'/%3e%3c/svg%3e");
}
.qc-app .form-switch .form-check-input:focus { box-shadow: 0 0 0 3px var(--qc-accent-soft); }
.qc-app .form-switch .form-check-label:empty,
.qc-app .form-switch .form-check-label > span:empty { display: none; }

/* -- Sliders (ion.rangeSlider, labels replaced by our own readout) -------- */
.qc-app .irs--shiny { height: 22px; }
.qc-app .irs--shiny .irs-line {
  top: 9px; height: 4px; border: 0; border-radius: 4px; background: var(--qc-track);
}
.qc-app .irs--shiny .irs-bar {
  top: 9px; height: 4px; border: 0; background: var(--qc-accent);
}
.qc-app .irs--shiny .irs-handle {
  top: 2px; width: 18px; height: 18px; border: 0; border-radius: 50%;
  background: #fff; cursor: grab;
  box-shadow: 0 0.5px 4px rgba(0, 0, 0, 0.35), 0 0 0 0.5px rgba(0, 0, 0, 0.06);
}
.qc-app .irs--shiny .irs-handle:hover, .qc-app .irs--shiny .irs-handle.state_hover { background: #fff; }
.qc-app .irs--shiny .irs-handle:active { cursor: grabbing; }
.qc-app .irs--shiny .irs-min, .qc-app .irs--shiny .irs-max,
.qc-app .irs--shiny .irs-from, .qc-app .irs--shiny .irs-to,
.qc-app .irs--shiny .irs-single, .qc-app .irs--shiny .irs-grid { display: none; }

/* -- Segmented controls --------------------------------------------------- */
/* Tabs (bslib navset_pill) */
.qc-main > .nav-pills, .qc-main .qc-tabs > .nav-pills, .qc-main > .tabbable > .nav-pills {
  display: flex; width: fit-content; flex-wrap: nowrap; gap: 2px; padding: 2px;
  border-radius: 10px; background: var(--qc-item); margin: 0 auto 16px;
}
.qc-main .nav-pills .nav-link {
  display: inline-flex; align-items: center; gap: 7px;
  padding: 6px 14px; border-radius: 8px;
  font-size: 12.5px; font-weight: 500; color: var(--qc-text-3);
  transition: background-color 0.15s ease, color 0.15s ease;
}
.qc-main .nav-pills .nav-link:hover { color: var(--qc-text); }
.qc-main .nav-pills .nav-link.active, .qc-main .nav-pills .show > .nav-link {
  background: var(--qc-seg-on); color: var(--qc-text);
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.04);
}
.qc-main .nav-pills .nav-link .fa, .qc-main .nav-pills .nav-link svg { opacity: 0.8; font-size: 11.5px; }

/* Radio-based segmented control (.qc_segmented) */
.qc-seg {
  display: inline-flex; gap: 2px; padding: 2px; margin: 0;
  border-radius: 9px; background: var(--qc-item); vertical-align: middle;
}
.qc-seg.shiny-input-container { width: auto !important; margin-bottom: 0; }
.qc-seg label { margin: 0; position: relative; display: inline-flex; }
.qc-seg input { position: absolute; opacity: 0; width: 1px; height: 1px; pointer-events: none; }
.qc-seg span {
  display: inline-flex; align-items: center; gap: 6px;
  padding: 4px 11px; border-radius: 7px; cursor: pointer;
  font-size: 12px; font-weight: 500; color: var(--qc-text-3); white-space: nowrap;
  transition: background-color 0.15s ease, color 0.15s ease;
}
.qc-seg span:hover { color: var(--qc-text); }
.qc-seg input:checked + span {
  background: var(--qc-seg-on); color: var(--qc-text);
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.04);
}
.qc-seg input:focus-visible + span { outline: 2px solid var(--qc-focus); outline-offset: 1px; }

/* Status summary in the tab row */
.qc-tab-status {
  position: absolute; top: 18px; right: 22px;
  display: flex; align-items: center; gap: 10px;
  font-size: 12px; color: var(--qc-text-3);
}
.qc-tab-status b { color: var(--qc-text); font-weight: 600; }

/* -- Cards ---------------------------------------------------------------- */
.qc-card {
  background: var(--qc-card); border: 1px solid var(--qc-border);
  border-radius: 16px; padding: 14px 16px; margin-bottom: 16px; min-width: 0;
}
.qc-card-head {
  display: flex; align-items: center; justify-content: space-between;
  gap: 12px; margin-bottom: 10px; min-height: 22px;
}
.qc-card-head-left { display: flex; align-items: baseline; gap: 10px; min-width: 0; flex-wrap: wrap; }
.qc-card-title {
  font-size: 11px; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.07em; color: var(--qc-text-3); white-space: nowrap;
}
.qc-card-sub { font-size: 12px; color: var(--qc-text-4); }
.qc-grid-2 { display: grid; grid-template-columns: minmax(0, 5fr) minmax(0, 7fr); gap: 16px; }
.qc-grid-2 > .qc-card { margin-bottom: 0; }
.qc-stack > * + * { margin-top: 16px; }

/* Toolbar rows (filters sit in one row above what they scope) */
.qc-toolbar {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap; margin-bottom: 12px;
}
.qc-toolbar .form-group { margin: 0; }
.qc-toolbar .qc-field { display: flex; align-items: center; gap: 7px; }
.qc-toolbar .qc-field > .qc-label { white-space: nowrap; }
.qc-toolbar .qc-spacer { flex: 1; }
.qc-hint { font-size: 11.5px; color: var(--qc-text-4); }
.qc-hint .fa { margin-right: 4px; }

/* Legend chips (HTML legends above charts) */
.qc-legend { display: flex; align-items: center; gap: 14px; font-size: 11.5px; color: var(--qc-text-3); flex-wrap: wrap; }
.qc-legend i { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; }
.qc-legend i.line { height: 0; width: 14px; border-radius: 0; border-top: 1.5px dashed var(--qc-text-2); vertical-align: 3px; }
.qc-legend i.pass { background: var(--qc-pass); }
.qc-legend i.fail { background: var(--qc-fail); }
.qc-legend i.fail-soft { background: var(--qc-fail-mid); }

/* -- Selects (selectize) -------------------------------------------------- */
.qc-app .selectize-control { margin: 0; }
.qc-app .selectize-input, .qc-app .selectize-input.full, .qc-app .selectize-input.focus {
  min-height: 30px; padding: 6px 28px 6px 10px !important;
  background: var(--qc-item); border: 0; border-radius: 8px; box-shadow: none;
  font-size: 12.5px; color: var(--qc-text); line-height: 18px;
}
.qc-app .selectize-input input { color: var(--qc-text); font-size: 12.5px; }
.qc-app .selectize-input input::placeholder { color: var(--qc-text-4); }
.qc-app .selectize-control.single .selectize-input::after { border-color: var(--qc-text-3) transparent transparent transparent; right: 11px; }
.qc-app .selectize-dropdown {
  background: var(--qc-menu); color: var(--qc-text-2);
  border: 1px solid var(--qc-border); border-radius: 10px;
  box-shadow: var(--qc-menu-shadow); margin-top: 4px; z-index: 1050; font-size: 12.5px;
}
.qc-app .selectize-dropdown-content { padding: 4px; max-height: 320px; }
.qc-app .selectize-dropdown .option, .qc-app .selectize-dropdown .optgroup-header {
  border-radius: 6px; padding: 6px 9px; color: var(--qc-text-2); background: transparent;
}
.qc-app .selectize-dropdown .optgroup-header {
  font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--qc-text-4);
}
.qc-app .selectize-dropdown .active { background: var(--qc-accent-soft); color: var(--qc-text); }
.qc-app .selectize-dropdown .selected { color: var(--qc-text); font-weight: 500; }
.qc-select-sm .selectize-control { min-width: 170px; }
.qc-select-cell .selectize-control { min-width: 240px; }
.qc-select-cell .selectize-input { font-family: var(--qc-mono); }

/* -- KPI tiles ------------------------------------------------------------ */
.qc-kpis { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-bottom: 16px; }
.qc-kpi {
  background: var(--qc-card); border: 1px solid var(--qc-border);
  border-radius: 16px; padding: 13px 16px 14px; min-width: 0;
}
.qc-kpi-label { display: flex; align-items: center; gap: 6px; font-size: 12px; color: var(--qc-text-3); }
.qc-kpi-label .fa { font-size: 11px; }
.qc-kpi-value { margin-top: 4px; font-size: 26px; font-weight: 600; letter-spacing: -0.02em; line-height: 1.15; color: var(--qc-text); }
.qc-kpi-value small { font-size: 13px; font-weight: 500; color: var(--qc-text-3); letter-spacing: 0; margin-left: 4px; }
.qc-kpi-sub { margin-top: 3px; font-size: 12px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-kpi-sub a { color: var(--qc-text-2); text-decoration: underline; text-underline-offset: 2px; cursor: pointer; }
.qc-kpi-sub a:hover { text-decoration: underline; }
.qc-kpi .icon-pass { color: var(--qc-good-text); }
.qc-kpi .icon-fail { color: var(--qc-bad-text); }
.qc-meter { height: 6px; border-radius: 3px; background: var(--qc-pass-soft); overflow: hidden; margin-top: 10px; }
.qc-meter > div { height: 100%; background: var(--qc-pass); border-radius: 3px; transition: width 0.3s ease; }

/* -- Pass / fail per metric -------------------------------------------- */
.qc-pf {
  display: grid; grid-template-columns: minmax(110px, max-content) minmax(60px, 1fr) auto auto;
  gap: 10px 16px; align-items: center;
}
.qc-pf-h {
  display: flex; align-items: center; gap: 6px;
  font-size: 10.5px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.06em;
  color: var(--qc-text-4); padding-bottom: 7px; border-bottom: 1px solid var(--qc-border);
}
.qc-pf-hn { justify-content: flex-end; }
.qc-pf-label { font-size: 12.5px; color: var(--qc-text-2); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
/* 100% bar: green share passes, red share fails, a 2 px gap between them */
.qc-pf-bar { display: flex; gap: 2px; height: 8px; }
.qc-pf-pass { background: var(--qc-pass); opacity: 0.55; border-radius: 4px 0 0 4px; }
.qc-pf-fail { background: var(--qc-fail); min-width: 3px; border-radius: 0 4px 4px 0; }
.qc-pf-bar > :only-child { border-radius: 4px; }
.qc-pf-n { font-size: 12.5px; color: var(--qc-text); text-align: right; min-width: 4.5em; font-variant-numeric: tabular-nums; }
.qc-pf-n.is-zero { color: var(--qc-text-4); }
/* the combined row, under a full-width rule */
.qc-pf-sep { grid-column: 1 / -1; border-top: 1px solid var(--qc-border); }
.qc-pf-label.is-total, .qc-pf-n.is-total { font-weight: 600; color: var(--qc-text); }
.qc-empty { padding: 26px 8px; text-align: center; font-size: 12.5px; color: var(--qc-text-4); }
.qc-empty .fa { display: block; font-size: 18px; margin-bottom: 8px; opacity: 0.7; }

/* -- Distributions grid --------------------------------------------------- */
.qc-dists { display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 10px; }
.qc-dist { background: var(--qc-inset); border-radius: 12px; padding: 9px 10px 4px; min-width: 0; }
.qc-dist-head { display: flex; justify-content: space-between; align-items: baseline; gap: 8px; }
.qc-dist-name { font-size: 12px; font-weight: 500; color: var(--qc-text-2); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-dist-gate { font-size: 11px; color: var(--qc-text-4); font-family: var(--qc-mono); white-space: nowrap; }
.qc-dist.is-off .qc-dist-name::after { content: "not used"; margin-left: 6px; font-size: 10.5px; color: var(--qc-text-4); font-weight: 400; }

/* -- Selection bar (lasso) ------------------------------------------------ */
.qc-selbar {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
  padding: 8px 10px 8px 14px; margin-bottom: 16px; border-radius: 12px;
  background: var(--qc-accent-soft); border: 1px solid var(--qc-border-2);
  font-size: 12.5px; color: var(--qc-text);
}
.qc-selbar .qc-spacer { flex: 1; }

/* -- Status pills --------------------------------------------------------- */
.qc-status {
  display: inline-flex; align-items: center; gap: 5px;
  padding: 3px 9px; border-radius: 999px; font-size: 11.5px; font-weight: 600; line-height: 1.2;
  white-space: nowrap;
}
.qc-status.pass { background: var(--qc-pass-soft); color: var(--qc-good-text); }
.qc-status.fail { background: var(--qc-fail-soft); color: var(--qc-bad-text); }
.qc-status.kept { background: transparent; color: var(--qc-good-text); box-shadow: inset 0 0 0 1px var(--qc-pass-line); }
.qc-status.excluded { background: transparent; color: var(--qc-bad-text); box-shadow: inset 0 0 0 1px var(--qc-fail-line); }
.qc-dot { display: inline-block; width: 7px; height: 7px; border-radius: 50%; flex: none; }
.qc-dot.pass { background: var(--qc-pass); }
.qc-dot.fail { background: var(--qc-fail); }
.qc-dot.off { background: var(--qc-text-4); }

/* -- Inspector ------------------------------------------------------------ */
.qc-cell-head { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.qc-cell-id { font-family: var(--qc-mono); font-size: 14px; font-weight: 600; color: var(--qc-text); }
.qc-cell-meta { font-size: 12px; color: var(--qc-text-3); }
.qc-cell-meta b { color: var(--qc-text-2); font-weight: 600; font-variant-numeric: tabular-nums; }
.qc-cell-fails { font-size: 12px; color: var(--qc-text-3); margin-top: 6px; }
.qc-cell-fails b { color: var(--qc-bad-text); font-weight: 500; }
.qc-pos { font-size: 12px; color: var(--qc-text-3); font-variant-numeric: tabular-nums; white-space: nowrap; }

.qc-tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(178px, 1fr)); gap: 10px; }
.qc-tile { background: var(--qc-inset); border-radius: 14px; padding: 10px 12px 9px; min-width: 0; }
.qc-tile-label { font-size: 11.5px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-tile-row { display: flex; align-items: baseline; justify-content: space-between; gap: 6px; }
.qc-tile-value { font-size: 18px; font-weight: 600; color: var(--qc-text); letter-spacing: -0.01em; }
.qc-tile-pct { font-size: 11px; color: var(--qc-text-4); font-variant-numeric: tabular-nums; white-space: nowrap; }
.qc-tile .qc-spark { height: 22px; margin: 6px 0 6px; }
.qc-tile-status { display: flex; align-items: center; gap: 6px; font-size: 11.5px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-tile-status.pass { color: var(--qc-good-text); }
.qc-tile-status.fail { color: var(--qc-bad-text); }
.qc-tile.is-fail { box-shadow: inset 0 0 0 1px var(--qc-fail-line); }

/* Copy-number legend above the profile plots */
.qc-cn-row { display: flex; justify-content: flex-end; margin: -2px 0 2px; }
.qc-cn-legend { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; font-size: 11px; color: var(--qc-text-3); }
.qc-cn-legend .qc-cn-title { color: var(--qc-text-4); }
.qc-cn-legend i { display: inline-block; width: 14px; height: 5px; border-radius: 2px; margin-right: 4px; vertical-align: 1px; }

.qc-explore-sel-head { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }

/* -- Tables (DT) ---------------------------------------------------------- */
.qc-app .dataTables_wrapper, .qc-app .dt-container { font-size: 12.5px; color: var(--qc-text-2); }
.qc-app table.dataTable { --bs-table-bg: transparent; --bs-table-color: var(--qc-text-2); --bs-table-striped-bg: transparent; border-collapse: collapse !important; }
.qc-app table.dataTable thead th {
  font-size: 10.5px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;
  color: var(--qc-text-3); white-space: nowrap; background: var(--qc-card);
  border-bottom: 1px solid var(--qc-border-2) !important; border-top: 0 !important;
}
.qc-app table.dataTable tbody td {
  border-top: 1px solid var(--qc-border) !important; border-bottom: 0 !important;
  font-variant-numeric: tabular-nums; white-space: nowrap; padding: 7px 10px;
  color: var(--qc-text-2); background: transparent; box-shadow: none !important;
}
.qc-app table.dataTable tbody td:first-child { font-family: var(--qc-mono); color: var(--qc-text); }
.qc-app table.dataTable tbody tr { cursor: pointer; }
.qc-app table.dataTable tbody tr:hover > * { background: var(--qc-hover) !important; }
.qc-app table.dataTable.no-footer { border-bottom: 0; }
.qc-app .dataTables_filter input, .qc-app .dt-search input,
.qc-app .dataTables_length select, .qc-app .dt-length select {
  background-color: var(--qc-item); border: 0; border-radius: 8px; color: var(--qc-text);
  font-size: 12.5px; padding: 5px 10px; box-shadow: none;
}
.qc-app .dataTables_info, .qc-app .dt-info { color: var(--qc-text-4); font-size: 12px; }
.qc-app .pagination { --bs-pagination-bg: transparent; --bs-pagination-border-color: transparent; --bs-pagination-color: var(--qc-text-3); --bs-pagination-hover-bg: var(--qc-item); --bs-pagination-hover-color: var(--qc-text); --bs-pagination-hover-border-color: transparent; --bs-pagination-active-bg: var(--qc-primary-bg); --bs-pagination-active-color: var(--qc-primary-text); --bs-pagination-active-border-color: var(--qc-primary-bg); --bs-pagination-disabled-bg: transparent; --bs-pagination-disabled-color: var(--qc-text-4); --bs-pagination-disabled-border-color: transparent; --bs-pagination-focus-box-shadow: none; gap: 2px; }
.qc-app .page-link { border-radius: 7px !important; font-size: 12px; padding: 4px 10px; }
.qc-app .dataTables_processing, .qc-app .dt-processing { background: var(--qc-menu); color: var(--qc-text-2); border-radius: 10px; border: 1px solid var(--qc-border); }

/* -- Tooltips, popovers, notifications, modals ---------------------------- */
.tooltip {
  --bs-tooltip-bg: var(--qc-menu); --bs-tooltip-color: var(--qc-text-2);
  --bs-tooltip-max-width: 270px; --bs-tooltip-font-size: 12px; --bs-tooltip-padding-x: 10px; --bs-tooltip-padding-y: 7px;
  --bs-tooltip-border-radius: 9px; --bs-tooltip-opacity: 1;
}
.tooltip .tooltip-inner { box-shadow: var(--qc-menu-shadow); text-align: left; line-height: 1.45; }
#shiny-notification-panel { width: 320px; }
.shiny-notification {
  background: var(--qc-menu); color: var(--qc-text-2); border: 1px solid var(--qc-border);
  border-radius: 12px; box-shadow: var(--qc-menu-shadow); font-size: 12.5px; opacity: 1;
}
.modal-content {
  --bs-modal-bg: var(--qc-card); --bs-modal-color: var(--qc-text); --bs-modal-border-color: var(--qc-border);
  --bs-modal-border-radius: 16px; --bs-modal-header-border-color: var(--qc-border); --bs-modal-footer-border-color: var(--qc-border);
}
.modal-title { font-size: 15px; font-weight: 600; }

/* -- Responsive ----------------------------------------------------------- */
@media (max-width: 1280px) {
  .qc-grid-2 { grid-template-columns: minmax(0, 1fr); }
}
@media (max-width: 1180px) {
  .qc-tab-status { position: static; justify-content: center; margin: 0 0 10px; }
}
@media (max-width: 1100px) {
  .qc-kpis { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .qc-pills .qc-pill:not(:first-child) { display: none; }
}
@media (max-width: 900px) {
  .qc-body, body.qc-sidebar-hidden .qc-body { grid-template-columns: minmax(0, 1fr); }
  .qc-sidebar { position: static; height: auto; border-right: 0; border-bottom: 1px solid var(--qc-border); }
  body.qc-sidebar-hidden .qc-sidebar { display: none; }
  .qc-main { padding: 14px 14px 40px; }
  .qc-title span, .qc-pills { display: none; }
}
@media (max-width: 560px) {
  .qc-kpis { grid-template-columns: minmax(0, 1fr); }
  .qc-main .nav-pills .nav-link { padding: 6px 9px; }
  .qc-main .nav-pills .nav-link .qc-tab-text { display: none; }
}

/* -- "Done" overlay (shown while the app hands the result back to R) ------ */
.qc-done {
  position: fixed; inset: 0; z-index: 3000;
  display: flex; align-items: center; justify-content: center; padding: 16px;
  background: var(--qc-glass);
  -webkit-backdrop-filter: blur(24px) saturate(180%);
  backdrop-filter: blur(24px) saturate(180%);
}
body:has(.qc-done) #shiny-disconnected-overlay { display: none; }
.qc-done-card {
  max-width: 420px; text-align: center; padding: 28px 30px 24px; border-radius: 20px;
  background: var(--qc-card); border: 1px solid var(--qc-border); box-shadow: var(--qc-menu-shadow);
}
.qc-done-card .qc-logo-hex { display: block; width: 104px; height: 120px; margin: 0 auto 14px; }
.qc-done-title { font-size: 17px; font-weight: 600; margin-bottom: 6px; }
.qc-done-text { font-size: 13px; color: var(--qc-text-3); line-height: 1.5; margin-bottom: 14px; }
.qc-done-card code {
  display: inline-block; font-family: var(--qc-mono); font-size: 12px; padding: 5px 10px;
  border-radius: 8px; background: var(--qc-item); color: var(--qc-text-2);
}

/* Shiny gives non-empty uiOutput() divs display: contents; these need a box. */
.qc-insp-head { display: block; width: 100%; margin-bottom: 10px; }
.qc-fails { color: var(--qc-text-2); }
.qc-more {
  display: inline-block; margin-left: 3px; padding: 1px 6px; border-radius: 999px;
  font-size: 10.5px; background: var(--qc-item); color: var(--qc-text-3);
}
)---"

# ── Client script (keyboard navigation, sidebar toggle) ──────────────────────
.qc_js <- r"---(
// ASCAT.sc QC app (launch_ui_qc): the only client-side behaviour.
// Everything else is driven from R.
(function () {
  // Left / right arrows (or k / j) step through cells in the Inspector.
  document.addEventListener('keydown', function (e) {
    if (e.metaKey || e.ctrlKey || e.altKey || e.defaultPrevented) return;
    var t = e.target;
    if (t && t.closest && t.closest('input, textarea, select, [contenteditable="true"], .irs, .selectize-control, .dropdown-menu')) return;
    var dir = 0;
    if (e.key === 'ArrowRight' || e.key === 'j') dir = 1;
    else if (e.key === 'ArrowLeft' || e.key === 'k') dir = -1;
    if (!dir || !window.Shiny || !Shiny.setInputValue) return;
    e.preventDefault();
    Shiny.setInputValue('qc_nav', dir, { priority: 'event' });
  });

  // Show / hide the gate sidebar; plotly charts re-measure afterwards.
  document.addEventListener('click', function (e) {
    var b = e.target.closest && e.target.closest('[data-qc-toggle-sidebar]');
    if (!b) return;
    var hidden = document.body.classList.toggle('qc-sidebar-hidden');
    b.setAttribute('aria-pressed', hidden ? 'true' : 'false');
    setTimeout(function () { window.dispatchEvent(new Event('resize')); }, 240);
  });
})();
)---"

# ── Logos: ASCAT.sc hex (circos v2), dark and light (lossless WebP, 240 px, base64)
.qc_logo_dark_b64 <- r"---(
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
txsPdbvxhLcdj3u78ewvzPV3ENjO36ngAA==
)---"

.qc_logo_light_b64 <- r"---(
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
qus/CIKc/1QwAQA=
)---"
