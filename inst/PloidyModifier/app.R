library(shiny)
library(ASCAT.sc)
library(shinybusy)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(shinyalert)
library(shinydashboard)
library(shinyFiles)
library(dipsaus)

reslocal <- NULL

pathrdata <- NULL
localOpt <- NULL

options(spinner.color="#f06313", spinner.color.background="#ffffff", spinner.size=1)

getIndex <- function(sample){
  # which() can return a vector if sample names are not unique;
  # [1] ensures we always get a scalar (NA when there is no match).
  index <- which(names(reslocal$allTracks.processed)==sample)
  return (index[1])
}

getSamples <- function() {
  if( !is.null(reslocal)){
    return (names(reslocal$allTracks.processed))}
}

#######################################################################
##################### QUALITY METRIC: MAPD ###############################
# MAPD (Median Absolute Pairwise Difference): the median of the absolute
# differences between consecutive bins' logR values on the UNSEGMENTED
# track -- a standard, assay-agnostic measure of bin-to-bin noise (lower
# is cleaner).
#
# Two data sources are needed depending on how the object was built:
#   - sc/shallow-coverage tracks (getTrackForAll()) keep the pre-
#     segmentation per-bin data in allTracks.processed[[sample]]$lSegs[[chr]]$data,
#     so compute_mapd() reads it straight from there.
#   - methylation-array tracks (getTrackForAll.bins()) discard that same
#     $data slot: $lSegs[[chr]]
#     only has $output (the segment table), $data is NULL, presumably to
#     save memory given an array carries hundreds of thousands of probes
#     vs. a few thousand bins for sc/shallow. For these,
#     compute_mapd_from_logr() reads the same pre-segmentation signal
#     from res$logr instead (the raw per-probe logR matrix, row-aligned
#     with res$annotations.probes) -- exactly what resegment_methyl_sample()
#     already reads successfully for the Re-segment tab.
# getAllSampleMAPD() tries the track-based method first and only falls
# back to the logr-based one if that comes up empty, so sc/shallow data
# is unaffected.
#######################################################################

compute_mapd <- function(track) {
  if (is.null(track) || is.null(track$lSegs)) return(NA_real_)

  all_diffs <- unlist(lapply(track$lSegs, function(seg) {
    d <- seg$data
    if (is.null(d) || nrow(d) < 2) return(NULL)
    if (!is.null(d$maploc)) d <- d[order(d$maploc), , drop = FALSE]

    # Identify the logR value column defensively: don't assume it's simply
    # "whatever's left after dropping chrom/maploc/chr/pos". Some track
    # builders carry an extra non-numeric column alongside the value column
    # (e.g. a probe/segment ID) -- blindly taking the first remaining
    # column can grab that instead, silently producing an all-NA vector.
    # Instead, try each remaining column in turn and use the first one
    # that's actually numeric.
    candidate_cols <- setdiff(colnames(d), c("chrom", "maploc", "chr", "pos"))
    vals <- NULL
    for (cc in candidate_cols) {
      v <- suppressWarnings(as.numeric(d[[cc]]))
      if (sum(is.finite(v)) >= 2) { vals <- v[is.finite(v)]; break }
    }
    if (is.null(vals) || length(vals) < 2) return(NULL)
    abs(diff(vals))
  }), use.names = FALSE)

  if (length(all_diffs) == 0) return(NA_real_)
  stats::median(all_diffs, na.rm = TRUE)
}

# Fallback for methylation objects: MAPD computed directly from the raw
# per-probe logR matrix (res$logr) rather than from allTracks.processed,
# since the latter's $data slot is discarded for methylation tracks. Note
# res$logr is the raw, PoN-normalised track -- *before* winsorising --
# so a handful of extreme unwinsorised probes could in principle appear;
# the median-based MAPD is robust to that (a few outlier diffs don't move
# a median), so this should still read out very close to what a
# winsorised-track MAPD would show.
compute_mapd_from_logr <- function(res, sample_index) {
  samp <- names(res$allTracks.processed)[sample_index]
  if (is.null(samp) || is.na(samp) || !samp %in% colnames(res$logr)) return(NA_real_)

  annot <- res$annotations.probes
  if (is.null(annot) || !all(c("chr", "pos") %in% colnames(annot))) return(NA_real_)

  vals_all  <- suppressWarnings(as.numeric(res$logr[, samp]))
  chrs      <- as.character(annot[, "chr"])
  positions <- suppressWarnings(as.numeric(as.character(annot[, "pos"])))

  if (length(vals_all) != length(chrs)) return(NA_real_)

  all_diffs <- unlist(lapply(split(seq_along(vals_all), chrs), function(idx) {
    ord <- idx[order(positions[idx])]
    v <- vals_all[ord]
    v <- v[is.finite(v)]
    if (length(v) < 2) return(NULL)
    abs(diff(v))
  }), use.names = FALSE)

  if (length(all_diffs) == 0) return(NA_real_)
  stats::median(all_diffs, na.rm = TRUE)
}

# MAPD for every currently loaded sample, named by sample. Cached on
# reslocal$mapd (a named list, sample name -> numeric MAPD or NA) so each
# sample is only ever computed once rather than on every Sample Viewer
# navigation, and so the values are preserved in any saved .Rda (they're
# a real field on the object, saved the same way as everything else).
# A sample's cached entry is explicitly dropped by the Re-segment/Refit
# "Apply" handlers whenever that sample's track is regenerated, so a
# stale MAPD can never be shown for a track that's since changed.
getAllSampleMAPD <- function() {
  samples <- getSamples()
  if (length(samples) == 0 || is.null(reslocal)) return(setNames(numeric(0), character(0)))

  cached  <- reslocal$mapd
  missing <- setdiff(samples, names(cached))

  if (length(missing) > 0) {
    has_raw_logr <- !is.null(reslocal$logr) && !is.null(reslocal$annotations.probes)
    new_vals <- lapply(missing, function(samp) {
      i  <- which(samples == samp)[1]
      tr <- tryCatch(reslocal$allTracks.processed[[i]], error = function(e) NULL)
      v  <- tryCatch(compute_mapd(tr), error = function(e) NA_real_)
      if (!is.na(v)) return(v)
      if (has_raw_logr) tryCatch(compute_mapd_from_logr(reslocal, i), error = function(e) NA_real_)
      else NA_real_
    })
    names(new_vals) <- missing
    reslocal$mapd <<- c(cached, new_vals)
    cached <- reslocal$mapd
  }

  vals <- unlist(cached[samples], use.names = TRUE)
  names(vals) <- samples
  vals
}

#######################################################################
##################### SAMPLE-VALIDITY GUARD ############################
# Centralises the "does this sample index actually have a usable      #
# solution / profile entry" check so every handler (not just View)    #
# can bail out cleanly instead of throwing on reslocal$X[[index]].    #
# Particularly important when the object holds a single sample: if    #
# that one sample has no valid entry, there is no other sample to     #
# fall back on, so every handler needs to fail gracefully.            #
#######################################################################

isSampleIndexValid <- function(index) {
  if (is.null(reslocal)) return(FALSE)
  if (is.null(index) || is.na(index)) return(FALSE)
  needed_lens <- c(
    length(reslocal$allTracks.processed),
    length(reslocal$allSolutions),
    length(reslocal$allSolutions.refitted.manual),
    length(reslocal$allProfiles.refitted.manual)
  )
  if (any(index > needed_lens)) return(FALSE)
  TRUE
}

warnInvalidSample <- function(sampleName, index) {
  shinyalert(
    "Warning",
    paste0("Sample '", sampleName, "' has no usable solution entry ",
           "(index ", index, "). Try a different sample, or reload the data."),
    type = "warning"
  )
}

#######################################################################
################ SAMPLE EXCLUSION / SUBSETTING ##########################
# Used by the Sample Viewer's "exclude sample" checkbox + "Save filtered
# .Rda" button. Simply doing res$allTracks.processed[keep_idx] is not
# enough: several other fields are indexed the exact same way and all
# have to be filtered in lockstep, or the saved object ends up with
# sample names mismatched against the wrong solution/profile. This
# covers every per-sample field the app itself reads/writes (including
# the manually-refitted profiles/solutions), plus, defensively, any
# other field that is clearly sample-keyed (a samples-as-columns
# matrix/data.frame such as a raw logR track, or a fully name-matched
# vector/list) even if this app never touches it directly.
#######################################################################

PER_SAMPLE_FIELDS <- c(
  "allTracks.processed", "allSolutions", "allProfiles",
  "allSolutions.refitted.auto",  "allProfiles.refitted.auto",
  "allSolutions.refitted.manual", "allProfiles.refitted.manual",
  "sex"
)

excludeSamplesFromRes <- function(res, exclude_names) {
  exclude_names <- unique(exclude_names)
  if (is.null(res) || length(exclude_names) == 0) return(res)

  all_samples <- names(res$allTracks.processed)
  if (is.null(all_samples) || length(all_samples) == 0) return(res)

  keep_names <- setdiff(all_samples, exclude_names)
  keep_idx   <- match(keep_names, all_samples)
  n_total    <- length(all_samples)

  # Filters one per-sample list/vector: match by name when possible
  # (safest -- correct even if the field is a different length, e.g.
  # `sex`), otherwise fall back to position only when the length lines
  # up exactly with the full original sample set.
  subsetOne <- function(val) {
    if (is.null(val)) return(val)
    if (!is.null(names(val)) && all(all_samples %in% names(val))) {
      val[keep_names]
    } else if (length(val) == n_total) {
      val[keep_idx]
    } else {
      val
    }
  }

  for (fld in PER_SAMPLE_FIELDS) {
    if (fld %in% names(res)) res[[fld]] <- subsetOne(res[[fld]])
  }

  # Anything else not already handled above: only touch it if it is
  # unambiguously sample-keyed, so unrelated fields (mode, gamma,
  # segmentation_alpha, etc.) are never accidentally altered.
  for (fld in setdiff(names(res), PER_SAMPLE_FIELDS)) {
    val <- res[[fld]]
    if (is.null(val)) next
    if (is.matrix(val) || is.data.frame(val)) {
      cn <- colnames(val)
      if (!is.null(cn) && any(cn %in% all_samples))
        res[[fld]] <- val[, setdiff(cn, exclude_names), drop = FALSE]
    } else if ((is.list(val) || is.atomic(val)) && !is.null(names(val)) &&
               all(all_samples %in% names(val))) {
      res[[fld]] <- val[keep_names]
    }
  }

  res
}

#######################################################################
######################### INTERFACE ###################################
#######################################################################

ui <- navbarPage(id="nav_page",
                 title="ASCAT.scFit", theme = shinytheme("united"),
                 tabPanel("Welcome",
                          busy_start_up(
                            loader = tags$img(src = "Loader.gif", width = 400),
                            text = "Loading ...",
                            color = "#ba4a00",
                            timeout = 3000,
                            background = "white",
                            mode = "auto"
                          ),
                          tags$head(
                            tags$style(HTML(
                              ".hover-effect {color: #FFFFFF ; border-radius: 30px; height:150px; width:300px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:4vh; border-width: 5px;  opacity: 0.90;}",
                              ".hover-effect:hover {color: #FFFFFF ; border-radius: 30px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:5vh; border-width: 5px; opacity: 1;}",
                              ".button-container {width: 300px;}",
                              ".click-effect:active {color: #FFFFFF ; border-radius: 30px; background-color: #3e0533; border-color: #6a0144; padding:5px; font-size:5vh; border-width: 5px; opacity: 1;}",
                              "code {display:block; padding:12px 16px; margin:auto; width:min(1100px, calc(100vw - 40px)); height:auto; max-height:min(500px, 72vh); overflow-y:auto; font-size:clamp(12px, 1.3vw, 14px); color:#3e0533; line-height:1.5; word-break:break-word; word-wrap:break-word; white-space:normal; background-color:#FFFFFF; border:6px solid #3e0533; border-radius:4px;}",
                              "pre {display:block; padding:9.5px; margin: auto; width: 350px; height: 100px; font-size:14px; color: #833e03; line-height:5px; word-break:break-all; word-wrap:break-word; white-space:pre-wrap; background-color:#fae6d4; border:4px solid #833e03; border-radius:4px;}",
                              ".sw-dropdown-content {  max-width: min(1200px, calc(100vw - 24px)) !important;  overflow-x: hidden;}",
                              # Shared button palette -- the same orange/purple used on this Welcome
                              # page (Load data's primary orange, Start's #3e0533/#6a0144 purple),
                              # reused across the Modifier / Re-segment / Refit / Sample Viewer tabs
                              # so every button in the app draws from the same two brand shades.
                              ".btn-app-orange {color:#FFFFFF !important; background-color:#ba4a00 !important; border-color:#ba4a00 !important; border-radius:8px; font-weight:bold; transition: background-color 0.15s ease, border-color 0.15s ease, transform 0.15s ease, box-shadow 0.15s ease;}",
                              ".btn-app-orange:hover {background-color:#e05a00 !important; border-color:#e05a00 !important; transform: translateY(-2px); box-shadow: 0 4px 10px rgba(186,74,0,0.35);}",
                              ".btn-app-orange:active {transform: translateY(0); box-shadow: none;}",
                              ".btn-app-purple {color:#FFFFFF !important; background-color:#3e0533 !important; border-color:#6a0144 !important; border-radius:8px; font-weight:bold; transition: background-color 0.15s ease, border-color 0.15s ease, transform 0.15s ease, box-shadow 0.15s ease;}",
                              ".btn-app-purple:hover {background-color:#6a0144 !important; border-color:#8a025c !important; transform: translateY(-2px); box-shadow: 0 4px 10px rgba(106,1,68,0.35);}",
                              ".btn-app-purple:active {transform: translateY(0); box-shadow: none;}"
                            ))),

                          #######################################################################
                          ######################### WELCOME TAB ################################
                          #######################################################################

                          img(src='WelcomeImage.png', align = "center",
                              style="width:100%; max-width:100%; position: absolute; z-index:-1;"),
                          div(style = "height:70px"),
                          fluidRow(
                            column(6, align="center",
                                   h1(strong("ASCAT.sc", style={'color: #5b016a; font-family: arial black ,sans-serif; font-size: 80px'}))),
                            column(6, align="center",
                                   h1(strong("Ploidy Modifier", style={'color: #5b016a; font-family: arial black ,sans-serif; font-size: 80px'})))
                          ),
                          br(), br(), br(), div(style = "height:50px"),

                          # ── Row 1: Help + compact Data type ──────────────────────────────
                          fluidRow(
                            column(6, align="center",
                                   dropdown(
                                     code(
                                       h4(strong("1. Choose the data", style={'font-family: Arial'}), align = "left"),
                                       p("Please choose the ASCAT.sc rdata object to load. Each sample in the object can be viewed and modified.", align = "left"),
                                       p("Please note that large datasets might take a while to load.", align = "left"),
                                       h4(strong("2. Select the data type", style={'font-family: Arial'}), align = "left"),
                                       p("Use the 'Data type' dropdown to specify whether your data comes from methylation arrays, single-cell sequencing, or shallow-coverage whole-genome sequencing.", align = "left"),
                                       h4(strong("3. Modify the profiles", style={'font-family: Arial; align-text: left'}), align = "left"),
                                       p("The purity & ploidy of the whole sample can be changed by clicking on the sunrise plot.", align = "left"),
                                       p("You can additionally shift the sample ploidy, or modify the copy number of 2 different segments.", align = "left"),
                                       p("Modifications can be done either by numeric input or directly on the Modifier Station plot.", align = "left"),
                                       h4(strong("4. Save the modified profiles", style={'font-family: Arial'}), align = "left"),
                                       p("Once you are happy with the results, you can save all the profiles in text format, or save the entire ASCAT.sc rdata object.", align = "left")
                                     ),
                                     size = "lg", circle = FALSE, status = "info",
                                     label = "Help", width = "1200px", inputId = "help"
                                   )
                            ),
                            # Data type: compact inline control styled to match the button row.
                            # A visible <select> drives a hidden Shiny selectInput so that
                            # input$data_type and updateSelectInput() work unchanged.
                            column(6, align="center",
                                   div(style = "display:inline-flex; align-items:center; gap:8px;
                                                background-color:rgba(255,255,255,0.80);
                                                border:2px solid #ba4a00;
                                                border-radius:6px; padding:5px 12px;",
                                       tags$span("Data type",
                                                 style = "font-family:Arial; font-size:13px;
                                                          color:#3e0533; font-weight:bold;
                                                          white-space:nowrap;"),
                                       # Visible, styled native <select>
                                       tags$select(
                                         id    = "data_type_vis",
                                         style = "border:1px solid #ccc; border-radius:4px;
                                                  padding:3px 6px; font-size:13px;
                                                  background:#fff; cursor:pointer;
                                                  height:30px; min-width:140px;",
                                         tags$option(value="methyl",  "Methylation"),
                                         tags$option(value="sc",      "Single Cell"),
                                         tags$option(value="shallow", "Shallow Coverage")
                                       ),
                                       # Hidden Shiny selectInput that keeps input$data_type in sync
                                       div(style = "display:none;",
                                           selectInput("data_type", label=NULL,
                                                       choices=c("Methylation"="methyl",
                                                                 "Single Cell"="sc",
                                                                 "Shallow Coverage"="shallow"),
                                                       selected="methyl", selectize=FALSE)
                                       ),
                                       # JS: visible -> hidden sync (both ways)
                                       tags$script(HTML(
                                         "$(document).on('change', '#data_type_vis', function(){
                                            $('#data_type').val($(this).val()).trigger('change');
                                          });
                                          $(document).on('shiny:inputchanged', function(e){
                                            if(e.name === 'data_type'){
                                              var v = e.value;
                                              if($('#data_type_vis').val() !== v)
                                                $('#data_type_vis').val(v);
                                            }
                                          });"
                                       ))
                                   )
                            )
                          ),

                          br(),

                          # ── Row 2: Load data – centre, primary CTA ───────────────────────
                          fluidRow(
                            column(4, align="center", offset = 4,
                                   shinyFilesButton("get_file", "Load data",
                                                    title    = "Choose the ASCAT.sc rdata object to load",
                                                    multiple = FALSE, buttonType = "primary",
                                                    class    = NULL,
                                                    style    = "font-size:16px; padding:10px 36px;
                                                                border-radius:8px; font-weight:bold;
                                                                width:200px;")
                            )
                          ),

                          br(), div(style = "height:30px"),

                          # ── Row 3: Start button ──────────────────────────────────────────────
                          fluidRow(
                            column(6, align="center", offset = 3,
                                   actionButtonStyled("start",
                                                      label  = "Start",
                                                      class  = "hover-effect click-effect",
                                                      width  = "300px")
                            )
                          )
                 ),

                 #######################################################################
                 ######################### MODIFIER TAB ################################
                 #######################################################################

                 tabPanel("Modifier",  shinyWidgets::useShinydashboard(),
                          fluidRow(box(width=12, title="Original", status="warning", solidHeader=TRUE,
                                       column(width=8, withSpinner(plotOutput("profile"), type=3)),
                                       column(width=4, plotOutput("sunrise1", height='450px')))
                          ),
                          useShinyjs(),
                          useShinyalert(),
                          fluidRow(
                            column(width=2, offset=1,
                                   dropdownButton(
                                     tags$h3("Choose Sample"),
                                     br(),
                                     selectInput(
                                       "samples", label = NULL,
                                       choices = getSamples(), selected = NULL,
                                       multiple = FALSE, selectize = FALSE
                                     ),
                                     size = "lg", circle = FALSE, status = "info",
                                     icon = icon("list", verify_fa = FALSE), width = "450px",
                                     label = "Select sample",
                                     actionButton("view", label = "View",
                                                  icon = icon("eye", verify_fa = FALSE),
                                                  class = "btn-app-orange")
                                   )
                            ),
                            column(width=2, offset=1,
                                   dropdownButton(
                                     tags$h3("Choose ploidy and purity on the sunrise graph"),
                                     br(),
                                     h5("To change the ploidy and purity of the whole sample directly on the Working station profile, please click the point on the sunrise plot corresponding to the desired ploidy/purity values. When 'Automatic optima' is enabled; the closest optimal purity/ploidy pair will be selected. Copy number segments can also be hidden to visualize the data points underneath. Neon red segments represent a copy number under 0 or above 8."),
                                     br(),
                                     checkboxInput("optima", "Automatic optima", TRUE),
                                     checkboxInput("hideCN", "Hide copy number segments", FALSE),
                                     checkboxInput("transparentCN", "Transparent copy number segments", FALSE),
                                     sliderInput("sizeP", "Increase point size:", min=1, max=3, value=1, step=0.5,
                                                 round=FALSE, ticks=TRUE, animate=FALSE, width=NULL,
                                                 sep=",", pre=NULL, post=NULL),
                                     size="lg", circle=FALSE, status="info",
                                     label="Modify purity & ploidy", width="450px"
                                   )
                            ),
                            column(width=5,
                                   dropdown(
                                     tags$h3("Choose a different way to modify the profile:"),
                                     dropdownButton(
                                       tags$h3("Modify the copy number of 2 segments"), br(),
                                       h5("The fit uses the LONGEST segment on each chosen chromosome."),
                                       selectInput("Chr1", "Choose first chromosome", c(1:22,"X","Y"),
                                                   selected=NULL, multiple=FALSE, selectize=FALSE, width="75%"),
                                       textInput("cn1", "Choose first copy number", value="", placeholder=NULL, width="75%"),
                                       selectInput("Chr2", "Choose second chromosome", c(1:22,"X","Y"),
                                                   selected=NULL, multiple=FALSE, selectize=FALSE, width="75%"),
                                       textInput("cn2", "Choose second copy number", value="", placeholder=NULL, width="75%"),
                                       actionButton("modify", label="Apply",
                                                    icon=icon("check", verify_fa=FALSE),
                                                    class="btn-app-orange"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Refit segments", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       tags$h3("Modify segment on graph"),
                                       h4("To modify the copy number of 2 different segments directly on the Working station profile, click anywhere on the desired chromosome, then click its target level. Repeat for a second (different) chromosome, then click 'Apply'. Note: the fit uses the LONGEST segment on each chromosome you click, not necessarily the exact segment under the cursor."),
                                       actionButton("refit", label="Apply",
                                                    icon=icon("check", verify_fa=FALSE),
                                                    class="btn-app-orange"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Refit segments on graph", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       sliderInput("ploidy", "Shift ploidy by:", -3, 4, 1, step=1,
                                                   round=FALSE, ticks=TRUE, animate=FALSE, width=NULL,
                                                   sep=",", pre=NULL, post=NULL),
                                       actionButton("shift", label="Apply",
                                                    icon=icon("check", verify_fa=FALSE),
                                                    class="btn-app-orange"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Shift sample ploidy", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       tags$h3("Shift ploidy on graph"),
                                       h4("To shift the ploidy of the whole sample directly on the Working station profile, please click on a point with the y axis position corresponding to the desired ploidy value, then click on 'Apply'"),
                                       actionButton("shift_graph", label="Apply",
                                                    icon=icon("check", verify_fa=FALSE),
                                                    class="btn-app-orange"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Shift ploidy on graph", width="450px", up=TRUE
                                     ),
                                     size="lg", circle=FALSE, status="info",
                                     label="Additional tools", width="800px", inputId="menu", up=TRUE
                                   )
                            )
                          ), br(),

                          fluidRow(
                            box(width=12, title="Modifier Working Station", status="warning", solidHeader=TRUE,
                                column(width=8, withSpinner(plotOutput("profile2", click="profile2_click"), type=3)),
                                column(width=4, plotOutput("sunrise2", click="sunrise2_click", height="460px"))
                            )
                          ),
                          fluidRow(
                            column(width=3, offset=1,
                                   actionButton("discard", label="Reset profile",
                                                icon=icon("arrows-rotate", verify_fa=FALSE),
                                                class="btn-app-orange")),
                            column(width=3, offset=1,
                                   downloadButton("savetxt", label="Save profiles",
                                                  class="btn-app-purple")),
                            column(width=3, offset=1,
                                   downloadButton("save", label="Save .Rda",
                                                  class="btn-app-purple"))
                          ), br(), br()
                 ),

                 #######################################################################
                 ######################### RE-SEGMENT TAB ###############################
                 #######################################################################
                 # Lets the user pick a sample and a different segmentation_alpha, and
                 # preview the profile that results from re-running everything after
                 # the raw logR extraction (winsorising, segmentation, purity/ploidy
                 # grid-search, profile fitting) with that new value -- without
                 # touching the Modifier working station until "Apply" is pressed.
                 #######################################################################

                 tabPanel("Re-segment", value="resegment_methyl_tab",
                          fluidRow(
                            column(width=12, align="left",
                                   h4(strong("Manual parameter adjustment",
                                             style={'font-family: Arial; color:#3e0533'})),
                                   p("Preview the effects of different parameters and if desired, save the profile.")
                            )
                          ),
                          fluidRow(
                            column(width=10, offset=1,
                              div(style="background-color:rgba(255,255,255,0.85); border:2px solid #ba4a00;
                                         border-radius:10px; padding:18px 22px; margin-bottom:16px;",
                                fluidRow(
                                  column(width=4,
                                         selectInput("resegment_sample", "Choose sample",
                                                     choices=getSamples(), selected=NULL,
                                                     multiple=FALSE, selectize=FALSE)
                                  ),
                                  column(width=4,
                                         selectInput("resegment_alpha", "Segmentation alpha",
                                                     choices = c("0.1","0.05","0.01","0.005",
                                                                 "0.001","0.0001","1e-05","1e-06",
                                                                 "1e-08","1e-10"),
                                                     selected = "0.001")
                                  ),
                                  column(width=4,
                                         sliderInput("resegment_gamma", "Gamma",
                                                     min = 0.4, max = 0.7, value = 0.55, step = 0.05,
                                                     ticks = TRUE)
                                  )
                                ),
                                fluidRow(
                                  column(width=7, style="padding-top:8px;",
                                         checkboxInput("resegment_quick",
                                                       "Quick preview",
                                                       value = TRUE)
                                  ),
                                  column(width=5, align="right",
                                         actionButton("resegment_run", label="Preview",
                                                      icon=icon("play", verify_fa=FALSE),
                                                      class="btn-app-orange")
                                  )
                                )
                              )
                            )
                          ),
                          fluidRow(
                            box(width=12, title="Re-segmented profile preview", status="warning", solidHeader=TRUE,
                                column(width=8, withSpinner(plotOutput("resegment_profile"), type=3)),
                                column(width=4, plotOutput("resegment_sunrise", height='450px')))
                          ),
                          fluidRow(
                            column(width=10, offset=1,
                              div(style="background-color:rgba(255,255,255,0.85); margin-bottom:16px;
                                         display:flex; align-items:center; gap:22px; flex-wrap:wrap;",
                                actionButton("resegment_apply",
                                             label="Apply to Modifier",
                                             icon=icon("check", verify_fa=FALSE),
                                             class="btn-app-purple"),
                                downloadButton("resegment_save", label="Save .Rda",
                                               class="btn-app-purple")
                              )
                            )
                          ), br(), br()
                 ),

                 #######################################################################
                 ######################### SC / SHALLOW REFIT TAB ########################
                 #######################################################################
                 # Same idea as the Re-segment tab above, but for single-cell / shallow-
                 # coverage data produced by run_sc_sequencing(): lets the user try a
                 # different segmentation_alpha AND a different binsize, re-binning the
                 # raw per-sample coverage track and re-segmenting it for one sample,
                 # without touching the Modifier working station until "Apply" is
                 # pressed.
                 #######################################################################

                 tabPanel("Refit (Single-Cell / Shallow)", value="resegment_sc_tab",
                          fluidRow(
                            column(width=12, align="left",
                                   h4(strong("Preview the effect of a different bin size and/or segmentation penalty",
                                             style={'font-family: Arial; color:#3e0533'})),
                                   p("Choose a sample, a new bin size, and a new value for segmentation_alpha (the DNAcopy segmentation penalty). This re-bins the raw per-sample coverage track and re-segments it for that sample only -- no re-reading of bams -- and previews the resulting profile below. Currently supported for single-cell / shallow-coverage sequencing data.")
                            )
                          ),
                          fluidRow(
                            column(width=10, offset=1,
                              div(style="background-color:rgba(255,255,255,0.85); border:2px solid #3e0533;
                                         border-radius:10px; padding:18px 22px; margin-bottom:16px;",
                                fluidRow(
                                  column(width=4,
                                         selectInput("screfit_sample", "Choose sample",
                                                     choices=getSamples(), selected=NULL,
                                                     multiple=FALSE, selectize=FALSE)
                                  ),
                                  column(width=4,
                                         numericInput("screfit_binsize", "Bin size (bp)",
                                                      value = 500000, min = 5000, step = 10000)
                                  ),
                                  column(width=4,
                                         selectInput("screfit_alpha", "Segmentation alpha",
                                                     choices = c("0.1","0.05","0.01","0.005",
                                                                 "0.001","0.0001","1e-05","1e-06",
                                                                 "1e-08","1e-10"),
                                                     selected = "0.01")
                                  )
                                ),
                                fluidRow(
                                  column(width=7, style="padding-top:8px;",
                                         checkboxInput("screfit_quick",
                                                       "Quick preview (coarser purity/ploidy grid, ~4x faster)",
                                                       value = TRUE)
                                  ),
                                  column(width=5, align="right",
                                         actionButton("screfit_run", label="Preview",
                                                      icon=icon("play", verify_fa=FALSE),
                                                      class="btn-app-orange")
                                  )
                                ),
                                fluidRow(
                                  column(width=12,
                                         textOutput("screfit_note"))
                                )
                              )
                            )
                          ),
                          fluidRow(
                            box(width=12, title="Refit preview", status="warning", solidHeader=TRUE,
                                column(width=8, withSpinner(plotOutput("screfit_profile"), type=3)),
                                column(width=4, plotOutput("screfit_sunrise", height='450px')))
                          ),
                          fluidRow(
                            column(width=10, offset=1,
                              div(style="background-color:rgba(255,255,255,0.85); border:2px solid #3e0533;
                                         border-radius:10px; padding:16px 22px; margin-bottom:16px;
                                         display:flex; align-items:center; gap:22px; flex-wrap:wrap;",
                                actionButton("screfit_apply",
                                             label="Apply to Modifier",
                                             icon=icon("check", verify_fa=FALSE),
                                             class="btn-app-purple"),
                                div(style="flex:1 1 320px; min-width:280px;",
                                    p(em("Applying updates the Modifier working station (allTracks.processed and allSolutions/allProfiles.refitted.manual) with the refitted result -- the original fit and any ASCAT.sc auto-refit are never overwritten. Note that Reset profile in the Modifier tab will still restore the original purity/ploidy values, but will display them against this newly re-binned/re-segmented track, since the segment boundaries themselves have changed."),
                                      style="color:#833e03; font-size:12px; line-height:1.5; margin:0;")
                                )
                              )
                            )
                          ), br(), br()
                 ),

                 #######################################################################
                 ######################### SAMPLE VIEWER TAB #############################
                 #######################################################################
                 # One profile at a time, full screen, with left/right navigation.
                 # No sunrise plot -- this page is a fast browsing tool, not an editor.
                 # Always shows the current "manual" (working) solution/profile for
                 # each sample, i.e. whatever is presently active in the Modifier tab.
                 #######################################################################

                 tabPanel("Sample Viewer",
                          tags$head(tags$style(HTML("
                            #sample_viewer_container { padding-top: 8px; }
                            #viewer_profile { height: 78vh; }
                            #viewer_profile img {
                              width: 100%; height: 100%;
                              object-fit: contain;
                            }
                            .viewer-arrow-btn {
                              font-size: 42px; line-height: 1; width: 100%; height: 78vh;
                              background: transparent; border: none; color: #ba4a00;
                              cursor: pointer; padding: 0;
                              transition: color 0.15s ease, transform 0.15s ease;
                            }
                            .viewer-arrow-btn:hover   { color: #6a0144; transform: scale(1.12); }
                            .viewer-arrow-btn:active  { transform: scale(0.94); }
                            .viewer-arrow-btn:disabled{ color: #dddddd; cursor: default; transform: none; }
                          "))),
                          div(id = "sample_viewer_container",
                              fluidRow(
                                column(width=6, offset=3, align="center",
                                       selectInput("viewer_sample_select", label=NULL,
                                                   choices=getSamples(), selected=NULL,
                                                   multiple=FALSE, selectize=TRUE, width="100%")
                                )
                              ),
                              fluidRow(
                                column(width=1, align="center",
                                       actionButton("viewer_prev_btn",
                                                    label=icon("chevron-left", verify_fa=FALSE),
                                                    class="viewer-arrow-btn")
                                ),
                                column(width=10, align="center",
                                       withSpinner(imageOutput("viewer_profile"), type=3)
                                ),
                                column(width=1, align="center",
                                       actionButton("viewer_next_btn",
                                                    label=icon("chevron-right", verify_fa=FALSE),
                                                    class="viewer-arrow-btn")
                                )
                              ),
                              fluidRow(
                                column(width=12, align="center",
                                       h5(textOutput("viewer_position_label"))
                                )
                              ),
                              fluidRow(
                                column(width=12, align="center",
                                       uiOutput("viewer_mapd")
                                )
                              ),
                              fluidRow(
                                column(width=12, align="center",
                                       tags$div(style="display:inline-block; padding:6px 14px; margin-top:4px;
                                                       background-color:rgba(255,255,255,0.85);
                                                       border:2px solid #ba4a00; border-radius:6px;",
                                                checkboxInput("viewer_exclude",
                                                              "Exclude this sample from saved output",
                                                              value = FALSE)
                                       )
                                )
                              ),
                              fluidRow(
                                column(width=10, offset=1,
                                  div(style="background-color:rgba(255,255,255,0.85); border:2px solid #ba4a00;
                                             border-radius:10px; padding:16px 22px; margin-top:14px;
                                             display:flex; align-items:center; justify-content:space-between;
                                             gap:18px; flex-wrap:wrap;",
                                    h5(textOutput("viewer_exclude_count"), style="color:#833e03; margin:0; white-space:nowrap;"),
                                    div(style="display:flex; align-items:center; gap:16px; flex-wrap:wrap;",
                                        actionButton("viewer_clear_exclusions", label="Clear exclusions",
                                                     icon=icon("broom", verify_fa=FALSE), class="btn-app-orange"),
                                        downloadButton("viewer_save", label="Save filtered .Rda", class="btn-app-purple")
                                    )
                                  )
                                )
                              ), br(),
                              # Left/right arrow-key navigation, active only while this tab
                              # is visible and the user isn't typing in a text/select field.
                              tags$script(HTML("
                                $(document).on('keydown', function(e) {
                                  if (!$('#sample_viewer_container').is(':visible')) return;
                                  var tag = (e.target.tagName || '').toLowerCase();
                                  if (tag === 'input' || tag === 'select' || tag === 'textarea') return;
                                  if (e.key === 'ArrowRight') {
                                    Shiny.setInputValue('viewer_next_key', Date.now(), {priority: 'event'});
                                    e.preventDefault();
                                  } else if (e.key === 'ArrowLeft') {
                                    Shiny.setInputValue('viewer_prev_key', Date.now(), {priority: 'event'});
                                    e.preventDefault();
                                  }
                                });
                              "))
                          )
                 )
)

#######################################################################
######################### SERVER ######################################
#######################################################################

server <- function(input, output, session) {

  volumes <- getVolumes()
  coords  <- reactiveValues(x=NULL, y=NULL)
  output$profile  <- renderPlot(NULL)
  output$profile2 <- renderPlot(NULL)
  output$resegment_profile <- renderPlot(NULL)
  output$resegment_sunrise <- renderPlot(NULL)
  output$screfit_profile   <- renderPlot(NULL)
  output$screfit_sunrise   <- renderPlot(NULL)

  # Only one of the two refit tabs is relevant for any given loaded
  # object (methylation vs sc/shallow); both start hidden and the
  # "start" observer below shows the right one once the mode is known.
  hideTab(inputId = "nav_page", target = "resegment_methyl_tab")
  hideTab(inputId = "nav_page", target = "resegment_sc_tab")

  # Keep only the last 4 clicks on the working-station profile: the
  # two-segment refit uses 4 (segment and target level, twice) and the
  # ploidy shift uses 1. The handlers clear this buffer as soon as they
  # read it.
  observeEvent(input$profile2_click, {
    coords$x <- utils::tail(c(coords$x, input$profile2_click$x), 4)
    coords$y <- utils::tail(c(coords$y, input$profile2_click$y), 4)
  })

  observeEvent(input$sunrise2_click, {
    coords$x <- input$sunrise2_click$x
    coords$y <- input$sunrise2_click$y
  })

  # Clicks belong to the profile they were made on: clear them when
  # another sample is selected.
  observeEvent(input$samples, {
    coords$x <- NULL
    coords$y <- NULL
  })

  sampleName <- reactive({ input$samples })
  optValue   <- reactive({ input$optima })
  shiftv     <- reactive({ input$ploidy })

  # ── Data type of the loaded object ──────────────────────────────────
  # Values: "sc" | "methyl" | "shallow". Decided once, when an object is
  # loaded (see detectDataMode() and the "start" observer), and stored in
  # loadedMode so that every handler reads the same value. Before anything
  # is loaded, the Welcome-page dropdown is used.
  loadedMode <- reactiveVal(NULL)
  dataMode <- reactive({
    if (!is.null(loadedMode())) loadedMode() else input$data_type
  })

  # Order of precedence: an explicit label on the object (res$mode,
  # res$is_methyl); then the object's structure (methylation results carry
  # the raw logR matrix and the probe annotations, sequencing results carry
  # the raw coverage tracks, the bin size or the genome build); then the
  # Welcome-page dropdown. For sequencing results without a label, the
  # dropdown only decides between single cell and shallow coverage.
  detectDataMode <- function(res, dropdown) {
    valid <- c("sc", "methyl", "shallow")
    if (length(res[["mode"]]) == 1 && res[["mode"]] %in% valid) return(res[["mode"]])
    if (isTRUE(res[["is_methyl"]])) return("methyl")
    if (!is.null(res[["logr"]]) && !is.null(res[["annotations.probes"]])) return("methyl")
    is_sequencing <- !is.null(res[["allTracks"]]) || !is.null(res[["binsize"]]) ||
      !is.null(res[["build"]])
    if (is_sequencing) return(if (identical(dropdown, "sc")) "sc" else "shallow")
    if (length(dropdown) == 1 && dropdown %in% valid) dropdown else "shallow"
  }

  # Euclidean distance between a click and a candidate optimum.
  point_dist <- function(x1, y1, x2, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2)

  cnHidden <- reactiveVal(FALSE)
  cnTrans  <- reactiveVal(FALSE)
  pSize    <- reactiveVal(FALSE)

  #########################################################################
  ######################## FILE CHOOSER ###################################
  #########################################################################

  observe({
    shinyFileChoose(input, "get_file", roots=volumes, session=session)
    if(!is.null(input$get_file)){
      file_selected <- parseFilePaths(volumes, input$get_file)
      pathrdata    <<- file_selected
    }
  })

  observe({ toggle(id="start", condition=(input$get_file >= 1)) })

  #########################################################################
  ######################## START ##########################################
  #########################################################################

  observeEvent(input$start, {
    updateNavbarPage(session=session, inputId="nav_page", selected="Sample Viewer")
    tryCatch({

      filepath <- as.character(pathrdata$datapath)
      load(filepath)
      reslocal <<- res

      # ── Defensive coercion for single-sample objects ─────────────────
      # If the object was built with sapply() over a single-sample loop
      # instead of lapply(), R may silently simplify a length-1 list into
      # a bare vector/object and/or drop its name. Everything downstream
      # (getSamples/getIndex/[[index]] lookups) assumes a *named list*,
      # so we normalise that here, right after loading, for every field
      # that is indexed per-sample elsewhere in the app.
      per_sample_fields <- c("allTracks.processed", "allSolutions", "allProfiles",
                              "allSolutions.refitted.auto", "allProfiles.refitted.auto",
                              "sex")
      for (fld in per_sample_fields) {
        if (fld %in% names(reslocal) && !is.null(reslocal[[fld]])) {
          if (!is.list(reslocal[[fld]]) || is.null(names(reslocal[[fld]]))) {
            if (!is.list(reslocal[[fld]]))
              reslocal[[fld]] <<- list(reslocal[[fld]])
            if (is.null(names(reslocal[[fld]])))
              names(reslocal[[fld]]) <<- if (length(reslocal[[fld]]) == 1)
                "Sample1" else seq_along(reslocal[[fld]])
          }
        }
      }

      updateSelectInput(session, "samples", label=NULL, choices=getSamples())
      updateSelectInput(session, "resegment_sample", label=NULL, choices=getSamples())
      updateSelectInput(session, "screfit_sample", label=NULL, choices=getSamples())
      updateSelectInput(session, "viewer_sample_select", choices=getSamples())
      viewerIndex(1L)
      excludedSamples(character(0))
      trackVersion(trackVersion() + 1L)

      # Work out the data type of the loaded object, store it for all
      # handlers, and show it in the Welcome-page dropdown.
      detected_mode <- detectDataMode(reslocal, isolate(input$data_type))
      loadedMode(detected_mode)
      updateSelectInput(session, "data_type", selected = detected_mode)

      mode <- detected_mode   # "sc" | "methyl" | "shallow"

      # ── Set gamma ──────────────────────────────────────────────────
      # Use the gamma stored on the loaded object when there is one
      # (res$gamma, or res$GAMMA); otherwise fall back to the default
      # for the data type.
      stored_gamma <- if (!is.null(reslocal[["gamma"]])) reslocal[["gamma"]] else reslocal[["GAMMA"]]
      stored_gamma <- suppressWarnings(as.numeric(stored_gamma))
      reslocal$gamma <<- if (length(stored_gamma) == 1 && is.finite(stored_gamma) && stored_gamma > 0)
        stored_gamma
      else if (mode == "methyl") 0.55 else 1

      # ── Show only the refit tab relevant to this data's mode ──────────
      if (mode == "methyl") {
        showTab(inputId = "nav_page", target = "resegment_methyl_tab")
        hideTab(inputId = "nav_page", target = "resegment_sc_tab")
        updateSliderInput(session, "resegment_gamma", value = reslocal$gamma)
      } else {
        hideTab(inputId = "nav_page", target = "resegment_methyl_tab")
        showTab(inputId = "nav_page", target = "resegment_sc_tab")
        updateNumericInput(session, "screfit_binsize",
                           value = if (!is.null(reslocal$binsize)) reslocal$binsize else 500000)
      }

      # ── Initialise manual slots if not already present ──────────────
      if (!"allSolutions.refitted.manual" %in% names(reslocal)) {
        if (mode == "sc") {
          reslocal$allProfiles.refitted.manual <<- reslocal$allProfiles
          reslocal$allSolutions.refitted.manual <<- reslocal$allSolutions
        } else {
          # methyl / shallow: prefer .refitted.auto if it exists
          if ("allSolutions.refitted.auto" %in% names(reslocal)) {
            reslocal$allProfiles.refitted.manual <<- reslocal$allProfiles.refitted.auto
            reslocal$allSolutions.refitted.manual <<- reslocal$allSolutions.refitted.auto
          } else {
            reslocal$allProfiles.refitted.manual <<- reslocal$allProfiles
            reslocal$allSolutions.refitted.manual <<- reslocal$allSolutions
          }
        }
      }

      # ── Make sure .refitted.auto exists for sc objects ──────────────
      # run_any_refitProfile() and run_any_refitProfile_shift() always
      # refit from res$allSolutions.refitted.auto[[sample]]. sc objects
      # only carry allSolutions/allProfiles, so seed the .refitted.auto
      # fields from those; the refit handlers below temporarily point
      # them at the current manual solution.
      if (mode == "sc") {
        if (!"allSolutions.refitted.auto" %in% names(reslocal))
          reslocal$allSolutions.refitted.auto <<- reslocal$allSolutions
        if (!"allProfiles.refitted.auto" %in% names(reslocal))
          reslocal$allProfiles.refitted.auto <<- reslocal$allProfiles
      }

      # ── Bail out early with a visible message if the (only) sample  ──
      # has no usable solution entry, rather than failing silently     #
      # further down or leaving blank plots with no explanation.       #
      if (!isSampleIndexValid(1)) {
        shinyalert("Error",
                   "The loaded data has no usable solution/profile entries for its samples. Please check the input object.",
                   type = "error")
        return(invisible(NULL))
      }

      # ── Original profile (top panel) ────────────────────────────────
      output$profile <- renderImage({
        outfile <- tempfile(fileext='.png')
        png(outfile, width=950, height=400)

        if (mode == "sc") {
          plotSolution(reslocal$allTracks.processed[[1]],
                       purity  = reslocal$allSolutions[[1]]$purity,
                       ploidy  = reslocal$allSolutions[[1]]$ploidy,
                       ismale  = resolve_ismale(reslocal$sex, 1),
                       gamma   = reslocal$gamma,
                       sol     = reslocal$allSolutions[[1]])
        } else {
          # methyl / shallow: use .refitted.auto if present
          if ("allSolutions.refitted.auto" %in% names(reslocal)) {
            plotSolution(reslocal$allTracks.processed[[1]],
                         purity = reslocal$allSolutions.refitted.auto[[1]]$purity,
                         ploidy = reslocal$allSolutions.refitted.auto[[1]]$ploidy,
                         ismale = resolve_ismale(reslocal$sex, 1),
                         gamma  = reslocal$gamma,
                         sol    = reslocal$allSolutions.refitted.auto[[1]])
          } else {
            plotSolution(reslocal$allTracks.processed[[1]],
                         purity = reslocal$allSolutions[[1]]$purity,
                         ploidy = reslocal$allSolutions[[1]]$ploidy,
                         ismale = resolve_ismale(reslocal$sex, 1),
                         gamma  = reslocal$gamma,
                         sol    = reslocal$allSolutions[[1]])
          }
        }

        dev.off()
        list(src=outfile, alt="Original profile", width='100%', height=400)
      }, deleteFile=TRUE)

      # ── Working-station profile (bottom panel) ───────────────────────
      output$profile2 <- renderPlot({
        isolate(plotSolution(reslocal$allTracks.processed[[1]],
                             purity        = reslocal$allSolutions.refitted.manual[[1]]$purity,
                             ploidy        = reslocal$allSolutions.refitted.manual[[1]]$ploidy,
                             ismale        = resolve_ismale(reslocal$sex, 1),
                             gamma         = reslocal$gamma,
                             sol           = reslocal$allSolutions.refitted.manual[[1]],
                             ambiguousFlag = FALSE,
                             hideCN        = cnHidden(),
                             transparentCN = cnTrans(),
                             zoomPoints    = pSize()))
      })

      # ── Sunrise plots ────────────────────────────────────────────────
      if (mode == "sc") {
        output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions[[1]], is_sc=TRUE)) })
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[1]], localMinima=TRUE, is_sc=TRUE)) })
      } else {
        if ("allSolutions.refitted.auto" %in% names(reslocal)) {
          output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions.refitted.auto[[1]])) })
        } else {
          output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions[[1]])) })
        }
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[1]], localMinima=TRUE)) })
      }

    }, error=function(e) {
      message("An Error Occurred while starting: ", conditionMessage(e))
      print(e)
      shinyalert("Error",
                 paste0("Could not load/initialise the data: ", conditionMessage(e)),
                 type = "error")
    })
  })

  #########################################################################
  ######################## DISCARD/KEEP ###################################
  #########################################################################

  observeEvent(input$discard, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
      mode  <- dataMode()

      if (mode == "sc") {
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity        = reslocal$allSolutions[[index]]$purity,
                               ploidy        = reslocal$allSolutions[[index]]$ploidy,
                               ismale        = resolve_ismale(reslocal$sex, index),
                               gamma         = reslocal$gamma,
                               sol           = reslocal$allSolutions[[index]],
                               ambiguousFlag = FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE)) })
        reslocal$allProfiles.refitted.manual[[index]]  <<- reslocal$allProfiles[[index]]
        reslocal$allSolutions.refitted.manual[[index]] <<- reslocal$allSolutions[[index]]

      } else {
        # methyl / shallow
        ref_sol  <- if ("allSolutions.refitted.auto" %in% names(reslocal)) reslocal$allSolutions.refitted.auto[[index]]  else reslocal$allSolutions[[index]]
        ref_prof <- if ("allProfiles.refitted.auto"  %in% names(reslocal)) reslocal$allProfiles.refitted.auto[[index]]   else reslocal$allProfiles[[index]]

        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity        = ref_sol$purity,
                               ploidy        = ref_sol$ploidy,
                               ismale        = resolve_ismale(reslocal$sex, index),
                               gamma         = reslocal$gamma,
                               sol           = ref_sol,
                               ambiguousFlag = FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE)) })
        reslocal$allProfiles.refitted.manual[[index]]  <<- ref_prof
        reslocal$allSolutions.refitted.manual[[index]] <<- ref_sol
      }
    } else {
      return(NULL)
    }
  })

  #########################################################################
  ######################## HIDE CN SEGMENTS ###############################
  #########################################################################

  observeEvent(input$hideCN, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
      cnHidden(input$hideCN)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = resolve_ismale(reslocal$sex, index),
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=input$hideCN, transparentCN=cnTrans(), zoomPoints=pSize()))
      })
    } else { return(NULL) }
  })

  observeEvent(input$transparentCN, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
      cnTrans(input$transparentCN)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = resolve_ismale(reslocal$sex, index),
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=cnHidden(), transparentCN=input$transparentCN, zoomPoints=pSize()))
      })
    } else { return(NULL) }
  })

  observeEvent(input$sizeP, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
      pSize(input$sizeP)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = resolve_ismale(reslocal$sex, index),
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=cnHidden(), transparentCN=input$transparentCN, zoomPoints=pSize()))
      })
    } else { return(NULL) }
  })

  #########################################################################
  ######################## DOWNLOAD #######################################
  #########################################################################

  output$save <- downloadHandler(
    filename = function() { "result_manualfitting.Rda" },
    content  = function(file) {
      showModal(modalDialog(div(tags$b("Loading...", style="color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      res <- reslocal
      save(res, file=file)
    }
  )

  # Same object/content as the Modifier tab's "Save .Rda" -- offered again
  # here so it's reachable without leaving the Re-segment page. Only
  # reflects whatever has already been "Applied" on this page (or
  # elsewhere); an un-applied preview isn't written into reslocal yet, so
  # it isn't included until Apply to Modifier is pressed first.
  output$resegment_save <- downloadHandler(
    filename = function() { "result_manualfitting.Rda" },
    content  = function(file) {
      showModal(modalDialog(div(tags$b("Loading...", style="color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      res <- reslocal
      save(res, file=file)
    }
  )

  output$savetxt <- downloadHandler(
    filename = function() { "Profiles_txt.zip" },
    content  = function(file) {
      showModal(modalDialog(div(tags$b("Loading...", style="color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      for (i in 1:length(reslocal$allProfiles.refitted.manual)) {
        write.table(reslocal$allProfiles.refitted.manual[[i]],
                    quote=F, sep="\t", col.names=T, row.names=F,
                    file=paste0(temp_directory, "/",
                                paste0(names(reslocal$allTracks.processed)[i], "-manual_refit"),
                                ".ASCAT.scprofile.txt"))
      }
      zip::zip(zipfile=file, files=dir(temp_directory), root=temp_directory)
    },
    contentType = "application/zip"
  )

  #########################################################################
  ######################## SUNRISE ########################################
  #########################################################################

  observeEvent(input$sunrise2_click, {
    index <- getIndex(sampleName())
    if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
    mode  <- dataMode()

    tryCatch({
      # Use the same solution (and error grid) that the working-station
      # sunrise plot displays, so clicks map onto what is drawn.
      solution <- reslocal$allSolutions.refitted.manual[[index]]
      errs     <- solution$errs
      errs     <- errs - min(errs)
      errs.max <- max(solution$errs[!is.infinite(solution$errs)])
      errs[is.infinite(errs)] <- errs.max
      errs <- errs / errs.max

      click_ploidy <- as.numeric(coords$x)
      click_purity <- as.numeric(coords$y)

      if (mode == "sc") {
        # SC uses a different row ordering convention
        purity_vals_click <- as.numeric(rownames(errs))
        if (purity_vals_click[1] < purity_vals_click[nrow(errs)]) {
          errs <- errs[rev(seq_len(nrow(errs))), ]
        }

        if (nrow(errs) == 1) {
          # ── Single-purity SC: only ploidy varies ──────────────────────
          if (optValue() && !is.null(localOpt) && !is.null(localOpt$bao)) {
            best <- which.min(abs(localOpt$bao[, 2] / ncol(errs) - click_ploidy))
            ploidy <- as.numeric(colnames(errs)[localOpt$bao[best, 2]])
          } else {
            col_idx <- max(1, min(ncol(errs), round(click_ploidy * ncol(errs))))
            ploidy  <- as.numeric(colnames(errs)[col_idx])
          }
          purity <- as.numeric(rownames(errs)[1])

        } else {
          # ── Normal multi-purity SC ────────────────────────────────────
          ploidy <- click_ploidy
          purity <- click_purity

          if (optValue() && !is.null(localOpt) && !is.null(localOpt$bao) &&
              nrow(localOpt$bao) > 0) {
            best <- 1
            dist <- point_dist(ploidy, purity,
                               localOpt$bao[best, 2] / ncol(errs),
                               1 - (localOpt$bao[best, 1] - 1) * 0.9 / (nrow(errs) - 1))
            for (i in seq_len(nrow(localOpt$bao))) {
              dist2 <- point_dist(ploidy, purity,
                                  localOpt$bao[i, 2] / ncol(errs),
                                  1 - (localOpt$bao[i, 1] - 1) * 0.9 / (nrow(errs) - 1))
              if (dist >= dist2) { dist <- dist2; best <- i }
            }
            ploidy <- localOpt$bao[best, 2] / ncol(errs)
            purity <- 1 - (localOpt$bao[best, 1] - 1) * 0.9 / (nrow(errs) - 1)
          }

          ploidy <- as.numeric(colnames(errs)[pmax(1L, pmin(ncol(errs), round(as.numeric(ploidy) * ncol(errs))))])
          purity <- as.numeric(rownames(errs)[
            max(1, min(nrow(errs), 1 + round((1 - as.numeric(purity)) * (nrow(errs) - 1) / 0.9, digits=0)))
          ])
        }

        reslocal$allProfiles.refitted.manual[[index]] <<- getProfile(
          fitProfile(tracksSingle=reslocal$allTracks.processed[[index]],
                     purity, ploidy, ismale=resolve_ismale(reslocal$sex, index), gamma=reslocal$gamma),
          CHRS=reslocal$chr)
        reslocal$allSolutions.refitted.manual[[index]]$ploidy <<- ploidy
        reslocal$allSolutions.refitted.manual[[index]]$purity <<- purity
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=purity, ploidy=ploidy, gamma=reslocal$gamma,
                               ismale=resolve_ismale(reslocal$sex, index),
                               sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
        })

      } else {
        # methyl / shallow share the same coordinate mapping. Reverse
        # rows only if needed -- mirrors plotSunrise()'s own (conditional)
        # orientation check exactly, so the click handler can never end up
        # reading a row order different from what's actually drawn.
        purity_vals_click <- as.numeric(rownames(errs))
        if (purity_vals_click[1] < purity_vals_click[nrow(errs)]) {
          errs <- errs[rev(seq_len(nrow(errs))), ]
        }

        if (nrow(errs) == 1) {
          # ── Single-purity methyl/shallow: only ploidy varies ──────────
          if (optValue() && !is.null(localOpt) && !is.null(localOpt$bao)) {
            best <- which.min(abs(localOpt$bao[, 2] / ncol(errs) - click_ploidy))
            ploidy <- as.numeric(colnames(errs)[localOpt$bao[best, 2]])
          } else {
            col_idx <- max(1, min(ncol(errs), round(click_ploidy * ncol(errs))))
            ploidy  <- as.numeric(colnames(errs)[col_idx])
          }
          purity <- as.numeric(rownames(errs)[1])

        } else {
          # ── Normal multi-purity methyl/shallow ────────────────────────
          ploidy <- click_ploidy
          purity <- click_purity

          if (optValue() && !is.null(localOpt) && !is.null(localOpt$bao) &&
              nrow(localOpt$bao) > 0) {
            best <- 1
            dist <- point_dist(ploidy, purity,
                               localOpt$bao[best, 2] / ncol(errs),
                               1 - localOpt$bao[best, 1] / nrow(errs))
            for (i in seq_len(nrow(localOpt$bao))) {
              dist2 <- point_dist(ploidy, purity,
                                  localOpt$bao[i, 2] / ncol(errs),
                                  1 - localOpt$bao[i, 1] / nrow(errs))
              if (dist >= dist2) { dist <- dist2; best <- i }
            }
            ploidy <- localOpt$bao[best, 2] / ncol(errs)
            purity <- 1 - localOpt$bao[best, 1] / nrow(errs)
          }

          ploidy <- as.numeric(colnames(errs)[pmax(1L, pmin(ncol(errs), round(as.numeric(ploidy) * ncol(errs))))])
          purity <- as.numeric(rownames(errs)[pmax(1L, pmin(nrow(errs), round((1 - as.numeric(purity)) * nrow(errs))))])
        }

        reslocal$allProfiles.refitted.manual[[index]] <<- getProfile(
          fitProfile(tracksSingle=reslocal$allTracks.processed[[index]],
                     purity, ploidy, ismale=resolve_ismale(reslocal$sex, index), gamma=reslocal$gamma),
          CHRS=reslocal$chr)
        reslocal$allSolutions.refitted.manual[[index]]$ploidy <<- ploidy
        reslocal$allSolutions.refitted.manual[[index]]$purity <<- purity
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=purity, ploidy=ploidy, gamma=reslocal$gamma,
                               ismale=resolve_ismale(reslocal$sex, index),
                               sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
        })
      }

      coords$y <- NULL
      coords$x <- NULL

    },
    error=function(e) {
      message('An Error Occurred'); print(e)
      shinyalert("Error", paste0("Cannot fit profile: ", conditionMessage(e)), type="error")
    })
  })

  #########################################################################
  ######################## MODIFY SEGMENTS ON GRAPH #######################
  #########################################################################

  observeEvent(input$refit, {
    index <- getIndex(sampleName())
    if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }

    # Four clicks are needed: segment 1 and its target level, then
    # segment 2 and its target level.
    if (length(coords$x) < 4 || length(coords$y) < 4) {
      shinyalert("Info",
                 "Please click on the segment then its target level, twice (once per segment), before pressing Apply.",
                 type = "info")
      return(NULL)
    }

    # Read and clear the click buffer so the points cannot be reused by
    # another click-driven tool.
    xs <- coords$x; ys <- coords$y
    coords$x <- NULL; coords$y <- NULL

    mode <- dataMode()
    # Chromosome boundaries must come from the SAME track that's actually
    # on screen (the currently selected sample), not always sample 1 --
    # different samples can have different segmentation/binning (e.g.
    # after an individual Re-segment) or a different chrY presence by sex,
    # any of which would silently shift every boundary after the mismatch.
    breaks <- c(0, cumsum(sapply(reslocal$allTracks.processed[[index]]$lSegs,
                                 function(x) max(x$output$loc.end)) / 1e+06))
    # Map the click's positional bin (1..length(breaks)-1) to the actual
    # chromosome label ("1".."22","X","Y") rather than passing the raw
    # position through as a number -- refitProfile() matches against
    # profile[,"chromosome"], which holds these labels, so a bare integer
    # like 23 never matches "X" and silently returns zero rows.
    chrLabels <- as.character(if (!is.null(reslocal$chr)) reslocal$chr else c(1:22, "X", "Y"))

    tryCatch({
      chr1 <- NULL; chr2 <- NULL
      for (i in seq_along(breaks)) { if (xs[1] < breaks[i]) { chr1 <- chrLabels[i-1]; break } }
      for (i in seq_along(breaks)) { if (xs[3] < breaks[i]) { chr2 <- chrLabels[i-1]; break } }
      if (is.null(chr1) || is.null(chr2) || is.na(chr1) || is.na(chr2)) {
        shinyalert("Error", "One of the clicks fell outside the plotted chromosomes -- please click again on the profile.", type="error")
        return(NULL)
      }
      if (identical(chr1, chr2)) {
        shinyalert("Info", "Both segments are on the same chromosome -- please click on two different chromosomes.", type="info")
        return(NULL)
      }
      y2 <- round(ys[2], digits=0)
      y4 <- round(ys[4], digits=0)

      if (mode == "sc") {
        rescopy <- reslocal
        # run_any_refitProfile() always refits FROM
        # res$allSolutions.refitted.auto[[sample_indice]], regardless of
        # mode -- so the field to point at the current manual solution
        # (and restore afterwards) must be .refitted.auto here too, not
        # bare allSolutions/allProfiles (which that function never reads).
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                        chr1=chr1, ind1=NA, n1=y2,
                                        chr2=chr2, ind2=NA, n2=y4,
                                        CHRS=reslocal$chr, outdir=tempdir(),
                                        gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
        })

      } else {
        # methyl / shallow
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                        chr1=chr1, ind1=NA, n1=y2,
                                        chr2=chr2, ind2=NA, n2=y4,
                                        CHRS=reslocal$chr, outdir=tempdir(),
                                        gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
        })
      }

      # refitProfile() applies the new solution even when the fit is
      # ambiguous, so report that flag to the user.
      if (isTRUE(reslocal$allSolutions.refitted.manual[[index]]$ambiguous)) {
        shinyalert("Warning",
                   "The new solution's fit was ambiguous; it was applied anyway -- please double check the result.",
                   type="warning")
      }

    },
    error=function(e) {
      message('An Error Occurred'); print(e)
      shinyalert("Error", paste0("Cannot fit profile: ", conditionMessage(e)), type="error")
    })
  })

  #########################################################################
  ######################## MODIFY SEGMENTS ################################
  #########################################################################

  observeEvent(input$modify, {
    index <- getIndex(sampleName())
    if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }

    # Chromosomes are passed as character labels ("1".."22","X","Y"),
    # matching profile[,"chromosome"] inside refitProfile().
    chr1v <- input$Chr1
    chr2v <- input$Chr2
    cn1v  <- suppressWarnings(as.numeric(input$cn1))
    cn2v  <- suppressWarnings(as.numeric(input$cn2))

    if (is.null(chr1v) || chr1v == "" || is.null(chr2v) || chr2v == "" ||
        is.na(cn1v) || is.na(cn2v)) {
      shinyalert("Info",
                 "Please choose both chromosomes and enter a numeric copy number for each segment.",
                 type = "info")
      return(NULL)
    }
    if (identical(as.character(chr1v), as.character(chr2v))) {
      shinyalert("Info",
                 "Please choose two different chromosomes: the fit needs one segment on each.",
                 type = "info")
      return(NULL)
    }

    mode <- dataMode()

    tryCatch({
      if (mode == "sc") {
        rescopy <- reslocal
        # See the "Modify segment on graph" handler above: run_any_refitProfile()
        # reads .refitted.auto regardless of mode, so that's what must be
        # swapped in/out here too (not bare allSolutions/allProfiles).
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                        chr1=chr1v, ind1=NA, n1=cn1v,
                                        chr2=chr2v, ind2=NA, n2=cn2v,
                                        CHRS=reslocal$chr, outdir=tempdir(),
                                        gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
        })

      } else {
        # methyl / shallow
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                        chr1=chr1v, ind1=NA, n1=cn1v,
                                        chr2=chr2v, ind2=NA, n2=cn2v,
                                        CHRS=reslocal$chr, outdir=tempdir(),
                                        gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
        })
      }

      if (isTRUE(reslocal$allSolutions.refitted.manual[[index]]$ambiguous)) {
        shinyalert("Warning",
                   "The new solution's fit was ambiguous; it was applied anyway -- please double check the result.",
                   type="warning")
      }

    },
    error=function(e) {
      message('An Error Occurred'); print(e)
      shinyalert("Error", paste0("Cannot fit profile: ", conditionMessage(e)), type="error")
    })
  })

  #########################################################################
  ######################## SHIFT ON GRAPH #################################
  #########################################################################

  observeEvent(input$shift_graph, {
    index <- getIndex(sampleName())
    if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }

    # At least one click on the profile is needed.
    if (length(coords$y) < 1) {
      shinyalert("Info", "Please click a point at the desired ploidy level on the profile before pressing Apply.", type = "info")
      return(NULL)
    }

    # Use the most recent click, and clear the buffer so it cannot be
    # reused.
    ypos <- as.numeric(round(utils::tail(coords$y, 1), digits=0))
    coords$x <- NULL; coords$y <- NULL

    mode <- dataMode()
    ploidy_cur <- round(as.numeric(reslocal$allSolutions.refitted.manual[[index]]$ploidy), digits=0)
    shiftp <- ypos - ploidy_cur

    tryCatch({
      if (mode == "sc") {
        rescopy <- reslocal
        # See "Modify segment on graph" above: run_any_refitProfile_shift()
        # always reads .refitted.auto, never bare allSolutions/allProfiles.
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                              CHRS=reslocal$chr, outdir=tempdir(),
                                              gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
        })

      } else {
        # methyl / shallow
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                              CHRS=reslocal$chr, outdir=tempdir(),
                                              gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                               ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                               ismale=resolve_ismale(reslocal$sex, index),
                               gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
        })
      }

      # refitProfile_shift() keeps the current solution and sets
      # $reverted when the new fit is ambiguous or not possible.
      if (isTRUE(reslocal$allSolutions.refitted.manual[[index]]$reverted)) {
        shinyalert("Warning", "New solution was ambiguous: reverted to the previous one.", type="warning")
      }

    },
    error=function(e) {
      message('An Error Occurred'); print(e)
      shinyalert("Error", paste0("Cannot fit profile: ", conditionMessage(e)), type="error")
    })
  })

  #########################################################################
  ######################## SHIFT ##########################################
  #########################################################################

  observeEvent(input$shift, {
    index <- getIndex(sampleName())
    if (!isSampleIndexValid(index)) { warnInvalidSample(sampleName(), index); return(NULL) }
    # This tool only reads the "Shift ploidy by" slider.
    mode   <- dataMode()
    shiftp <- as.numeric(shiftv())

    tryCatch({
      if (mode == "sc") {
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                              CHRS=reslocal$chr, outdir=tempdir(),
                                              gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          plotSolution(reslocal$allTracks.processed[[index]],
                       purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                       ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                       ismale=resolve_ismale(reslocal$sex, index),
                       gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                       ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize())
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
        })

      } else {
        # methyl / shallow
        rescopy <- reslocal
        rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
        rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
        rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                              CHRS=reslocal$chr, outdir=tempdir(),
                                              gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
        rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
        rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
        reslocal <<- rescopy
        output$profile2 <- renderPlot({
          plotSolution(reslocal$allTracks.processed[[index]],
                       purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                       ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                       ismale=resolve_ismale(reslocal$sex, index),
                       gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                       ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize())
        })
        output$sunrise2 <- renderPlot({
          localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE)
        })
      }

      if (isTRUE(reslocal$allSolutions.refitted.manual[[index]]$reverted)) {
        shinyalert("Warning", "New solution was ambiguous: reverted to the previous one.", type="warning")
      }

    },
    error=function(e) {
      shinyalert("Error", paste0("Cannot fit profile: ", conditionMessage(e)), type="error")
    })
  })

  #########################################################################
  ######################## VIEW ###########################################
  #########################################################################

  observeEvent(input$view, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      mode  <- dataMode()

      # Guard: allSolutions may be shorter than allTracks.processed if some
      # samples had no valid solution in the original run. Uses the shared
      # isSampleIndexValid() helper so the same check applies everywhere.
      if (!isSampleIndexValid(index)) {
        warnInvalidSample(sampleName(), index)
        return(NULL)
      }

      output$profile <- renderImage({
        outfile <- tempfile(fileext='.png')
        png(outfile, width=950, height=400)

        if (mode == "sc") {
          plotSolution(reslocal$allTracks.processed[[index]],
                       purity=reslocal$allSolutions[[index]]$purity,
                       ploidy=reslocal$allSolutions[[index]]$ploidy,
                       ismale=resolve_ismale(reslocal$sex, index),
                       gamma=reslocal$gamma, sol=reslocal$allSolutions[[index]])
        } else {
          # methyl / shallow
          if ("allSolutions.refitted.auto" %in% names(reslocal)) {
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions.refitted.auto[[index]]$purity,
                         ploidy=reslocal$allSolutions.refitted.auto[[index]]$ploidy,
                         ismale=resolve_ismale(reslocal$sex, index),
                         gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.auto[[index]])
          } else {
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions[[index]]$purity,
                         ploidy=reslocal$allSolutions[[index]]$ploidy,
                         ismale=resolve_ismale(reslocal$sex, index),
                         gamma=reslocal$gamma, sol=reslocal$allSolutions[[index]])
          }
        }

        dev.off()
        list(src=outfile, alt="Original profile", width='100%', height=400)
      }, deleteFile=TRUE)

      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = resolve_ismale(reslocal$sex, index),
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
      })

      if (mode == "sc") {
        output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions[[index]], is_sc=TRUE)) })
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE)) })
      } else {
        # methyl / shallow
        if ("allSolutions.refitted.auto" %in% names(reslocal)) {
          output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions.refitted.auto[[index]])) })
        } else {
          output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions[[index]])) })
        }
        output$sunrise2 <- renderPlot({ isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE)) })
      }

    } else { return(NULL) }
  })

  #########################################################################
  ######################## RE-SEGMENT #####################################
  #########################################################################
  # Preview: rerun winsorising -> segmentation -> purity/ploidy fitting ->
  # profile for ONE sample, using a new segmentation_alpha, starting from
  # the raw logR track already stored in reslocal$logr. This does NOT
  # touch reslocal until "Apply to Modifier" is pressed.
  #########################################################################

  resegmentPreview <- reactiveVal(NULL)

  # Shows the sample's CURRENT fit (the same "auto"-preferred-over-original
  # convention used by the View button) in the preview panel, so it's never
  # blank before "Preview" is run -- Preview then overwrites both plots with
  # the new segmentation_alpha/gamma result, same as before.
  showResegmentBaseline <- function(sample_name) {
    index <- getIndex(sample_name)
    if (!isSampleIndexValid(index)) return(invisible(NULL))

    base_sol <- if ("allSolutions.refitted.auto" %in% names(reslocal))
                  reslocal$allSolutions.refitted.auto[[index]]
                else
                  reslocal$allSolutions[[index]]
    if (is.null(base_sol)) return(invisible(NULL))

    ismale <- resolve_ismale(reslocal$sex, index, sample_name)

    output$resegment_profile <- renderPlot({
      isolate(plotSolution(reslocal$allTracks.processed[[index]],
                           purity = base_sol$purity,
                           ploidy = base_sol$ploidy,
                           ismale = ismale,
                           gamma  = reslocal$gamma,
                           sol    = base_sol))
    })
    output$resegment_sunrise <- renderPlot({
      isolate(plotSunrise(base_sol))
    })
  }

  # `reslocal` is a plain global (not a reactiveVal), so reading it inside
  # observe() creates no dependency and the block would only ever run once,
  # at app start-up, before any data is loaded. Refresh explicitly instead:
  # once right after loading (see updateSelectInput call in the "start"
  # observer above), and again defensively every time this tab is opened.
  observeEvent(input$nav_page, {
    if (identical(input$nav_page, "resegment_methyl_tab")) {
      updateSelectInput(session, "resegment_sample", choices = getSamples())
      if (!is.null(input$resegment_sample) && input$resegment_sample != "")
        showResegmentBaseline(input$resegment_sample)
    }
  })

  observeEvent(input$resegment_sample, {
    req(input$resegment_sample)
    showResegmentBaseline(input$resegment_sample)
  }, ignoreInit = TRUE)

  observeEvent(input$resegment_run, {
    req(input$resegment_sample)

    if (is.null(reslocal)) {
      shinyalert("Error", "Please load data first.", type="error")
      return(NULL)
    }

    # Gate on the fields the resegmentation actually needs, rather than on
    # a "mode"/"is_methyl" label: run_methylation_array()'s output can pass
    # through predictRefit_all()/printResults_all() (ASCAT.sc internals we
    # don't control), which may or may not preserve custom top-level fields
    # like is_methyl when they rebuild the result list. Checking directly
    # for res$logr and res$annotations.probes is what the function itself
    # requires, so it can't go stale the way a label-based check can, and
    # it tells us exactly what's missing if it fails.
    required_fields <- c("logr", "annotations.probes")
    missing_fields <- required_fields[
      !required_fields %in% names(reslocal) |
      vapply(required_fields, function(f) is.null(reslocal[[f]]), logical(1))
    ]
    if (length(missing_fields) > 0) {
      shinyalert("Not supported for this data",
                 paste0("Segmentation-alpha preview needs res$logr and res$annotations.probes ",
                        "(the raw logR track and probe annotations produced by run_methylation_array()). ",
                        "This loaded object is missing: ", paste(missing_fields, collapse=", "), ". ",
                        "This feature currently only supports methylation-array results that retain those fields."),
                 type="warning")
      return(NULL)
    }

    index <- getIndex(input$resegment_sample)
    if (!isSampleIndexValid(index)) { warnInvalidSample(input$resegment_sample, index); return(NULL) }

    alpha <- as.numeric(input$resegment_alpha)
    if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
      shinyalert("Invalid value", "segmentation_alpha must be a number strictly between 0 and 1.", type="warning")
      return(NULL)
    }

    gamma_choice <- input$resegment_gamma
    if (is.na(gamma_choice) || gamma_choice < 0.4 || gamma_choice > 0.7) {
      shinyalert("Invalid value", "gamma must be between 0.4 and 0.7.", type="warning")
      return(NULL)
    }

    showModal(modalDialog(div(tags$b("Re-segmenting and refitting...", style="color: steelblue;")), footer=NULL))

    # Quick preview trades grid resolution for speed: the full default grid
    # (purs step 0.01 x ploidies step 0.01) is ~21,000 points per sample;
    # doubling the step size cuts that to ~1/4 with a negligible effect on
    # which purity/ploidy region gets picked for a *preview*.
    grid_step <- if (isTRUE(input$resegment_quick)) 0.02 else 0.01
    purs_grid     <- seq(0.1, 1,   grid_step)
    ploidies_grid <- seq(1.7, 4,   grid_step)

    resegment_warnings <- character(0)
    out <- tryCatch(
      withCallingHandlers(
        resegment_methyl_sample(reslocal, sample_index = index, segmentation_alpha = alpha,
                                 gamma = gamma_choice, purs = purs_grid, ploidies = ploidies_grid),
        warning = function(w) {
          resegment_warnings <<- c(resegment_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) { message("Resegment error: ", conditionMessage(e)); print(e); NULL }
    )
    if (length(resegment_warnings) > 0)
      message("Resegment warnings (", length(resegment_warnings), "): ",
              paste(unique(resegment_warnings), collapse=" | "))

    removeModal()

    if (is.null(out)) {
      shinyalert("Error", "Could not re-segment/re-fit this sample with the chosen alpha. Please choose a different value.", type="error")
      return(NULL)
    }

    resegmentPreview(list(index = index, sample = input$resegment_sample,
                          params = out$params, result = out))

    ismale <- resolve_ismale(reslocal$sex, index, input$resegment_sample)

    output$resegment_profile <- renderPlot({
      isolate(plotSolution(out$track,
                           purity = out$solution$purity,
                           ploidy = out$solution$ploidy,
                           ismale = ismale,
                           gamma  = out$params$gamma,
                           sol    = out$solution))
    })
    output$resegment_sunrise <- renderPlot({
      isolate(plotSunrise(out$solution))
    })
  })

  observeEvent(input$resegment_apply, {
    preview <- resegmentPreview()
    if (is.null(preview)) {
      shinyalert("Nothing to apply", "Run a preview first.", type="warning")
      return(NULL)
    }

    index <- preview$index
    out   <- preview$result

    # Attach the exact parameters this fit used directly onto the solution,
    # so a record of them travels with the profile (visible on inspection,
    # and preserved in any subsequently saved .Rda/.txt).
    stamped_solution <- annotate_solution_with_params(out$solution, preview$params)

    rescopy <- reslocal
    # The re-segmented TRACK has to replace allTracks.processed -- the new
    # solution/profile were fit against it, so leaving the old track in
    # place would show new purity/ploidy values over stale segment
    # boundaries. The ORIGINAL fit (allSolutions/allProfiles) and the
    # ASCAT.sc "auto" refit (allSolutions.refitted.auto/
    # allProfiles.refitted.auto) are intentionally left untouched -- only
    # the "manual" working-station fields are ever written here, so
    # Reset/Discard in the Modifier tab can still recover the pre-refit
    # baseline for the purity/ploidy solution.
    rescopy$allTracks.processed[[index]] <- out$track
    rescopy$allSolutions.refitted.manual[[index]] <- stamped_solution
    rescopy$allProfiles.refitted.manual[[index]]  <- out$profile
    rescopy$refit_log <- append_refit_log(rescopy$refit_log, preview$sample, "methyl", preview$params)
    # The track just changed, so any cached MAPD for this sample is stale --
    # drop it (not recompute here) so it's simply recomputed, once, next
    # time it's actually needed.
    if (!is.null(rescopy$mapd) && preview$sample %in% names(rescopy$mapd))
      rescopy$mapd[[preview$sample]] <- NULL
    reslocal <<- rescopy
    trackVersion(trackVersion() + 1L)

    shinyalert("Applied",
               paste0("The re-segmented profile for '", preview$sample,
                      "' (segmentation_alpha=", preview$params$segmentation_alpha,
                      ", gamma=", preview$params$gamma,
                      ") has been applied to the Modifier working station (allSolutions.refitted.manual / ",
                      "allProfiles.refitted.manual). The original fit is untouched. ",
                      "Open the Modifier tab to continue refining it."),
               type="success")
  })

  #########################################################################
  ######################## SC / SHALLOW REFIT ##############################
  #########################################################################
  # Preview: re-bin the raw per-sample coverage track at a new binsize,
  # then re-segment/re-fit with a new segmentation_alpha, for ONE sample
  # of an sc/shallow-coverage result (run_sc_sequencing() output). This
  # does NOT touch reslocal until "Apply to Modifier" is pressed.
  #########################################################################

  screfitPreview <- reactiveVal(NULL)

  observeEvent(input$nav_page, {
    if (identical(input$nav_page, "resegment_sc_tab"))
      updateSelectInput(session, "screfit_sample", choices = getSamples())
  })

  output$screfit_note <- renderText({ "" })

  observeEvent(input$screfit_run, {
    req(input$screfit_sample)

    if (is.null(reslocal)) {
      shinyalert("Error", "Please load data first.", type="error")
      return(NULL)
    }

    # Gate on the fields the re-binning/re-segmentation actually needs,
    # same rationale as the methylation Re-segment tab: check directly
    # for what resegment_sc_sample() requires rather than trusting a
    # "mode" label.
    required_fields <- c("allTracks", "build", "chr")
    missing_fields <- required_fields[
      !required_fields %in% names(reslocal) |
      vapply(required_fields, function(f) is.null(reslocal[[f]]), logical(1))
    ]
    if (length(missing_fields) > 0) {
      shinyalert("Not supported for this data",
                 paste0("Bin size / segmentation-alpha preview needs res$allTracks (raw per-sample ",
                        "coverage), res$build, and res$chr, as produced by run_sc_sequencing(). ",
                        "This loaded object is missing: ", paste(missing_fields, collapse=", "), "."),
                 type="warning")
      return(NULL)
    }

    index <- getIndex(input$screfit_sample)
    if (!isSampleIndexValid(index)) { warnInvalidSample(input$screfit_sample, index); return(NULL) }

    samp <- names(reslocal$allTracks.processed)[index]
    if (is.null(reslocal$allTracks[[samp]]) || is.null(reslocal$allTracks[[samp]]$lCTS.tumour)) {
      shinyalert("Not supported for this sample",
                 paste0("res$allTracks[['", samp, "']]$lCTS.tumour (the raw, unbinned coverage track) ",
                        "was not found -- cannot re-bin this sample."),
                 type="warning")
      return(NULL)
    }

    alpha <- as.numeric(input$screfit_alpha)
    if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
      shinyalert("Invalid value", "segmentation_alpha must be a number strictly between 0 and 1.", type="warning")
      return(NULL)
    }

    binsize <- input$screfit_binsize
    if (is.na(binsize) || binsize <= 0) {
      shinyalert("Invalid value", "binsize must be a positive number.", type="warning")
      return(NULL)
    }

    showModal(modalDialog(div(tags$b("Re-binning, re-segmenting and refitting...", style="color: steelblue;")), footer=NULL))

    # Quick preview trades grid resolution for speed, same rationale as
    # the methylation tab.
    grid_step <- if (isTRUE(input$screfit_quick)) 0.02 else 0.01
    purs_grid     <- seq(0.1, 1, grid_step)
    ploidies_grid <- seq(1.7, 5, grid_step)

    screfit_warnings <- character(0)
    out <- tryCatch(
      withCallingHandlers(
        resegment_sc_sample(reslocal, sample_index = index, segmentation_alpha = alpha,
                             binsize = binsize, purs = purs_grid, ploidies = ploidies_grid),
        warning = function(w) {
          screfit_warnings <<- c(screfit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) { message("SC refit error: ", conditionMessage(e)); print(e); NULL }
    )
    if (length(screfit_warnings) > 0)
      message("SC refit warnings (", length(screfit_warnings), "): ",
              paste(unique(screfit_warnings), collapse=" | "))

    removeModal()

    if (is.null(out)) {
      shinyalert("Error", "Could not re-bin/re-segment/re-fit this sample with the chosen settings. Please choose different values.", type="error")
      return(NULL)
    }

    screfitPreview(list(index = index, sample = input$screfit_sample,
                        params = out$params, result = out))

    ismale  <- resolve_ismale(reslocal$sex, index, input$screfit_sample)
    is_sc_m <- identical(dataMode(), "sc")

    output$screfit_profile <- renderPlot({
      isolate(plotSolution(out$track,
                           purity = out$solution$purity,
                           ploidy = out$solution$ploidy,
                           ismale = ismale,
                           gamma  = reslocal$gamma,
                           sol    = out$solution))
    })
    output$screfit_sunrise <- renderPlot({
      isolate(plotSunrise(out$solution, is_sc = is_sc_m))
    })
    output$screfit_note <- renderText({ if (!is.null(out$note)) out$note else "" })
  })

  observeEvent(input$screfit_apply, {
    preview <- screfitPreview()
    if (is.null(preview)) {
      shinyalert("Nothing to apply", "Run a preview first.", type="warning")
      return(NULL)
    }

    index <- preview$index
    out   <- preview$result

    # See the equivalent comment in the methylation resegment_apply
    # handler: the track has to be replaced (the new solution was fit
    # against it), but the original fit and the ASCAT.sc "auto" refit are
    # never touched -- only the "manual" working-station fields are
    # written here.
    stamped_solution <- annotate_solution_with_params(out$solution, preview$params)

    rescopy <- reslocal
    rescopy$allTracks.processed[[index]] <- out$track
    rescopy$allSolutions.refitted.manual[[index]] <- stamped_solution
    rescopy$allProfiles.refitted.manual[[index]]  <- out$profile
    rescopy$refit_log <- append_refit_log(rescopy$refit_log, preview$sample, "sc", preview$params)
    # The track just changed (new binsize and/or segmentation_alpha), so
    # any cached MAPD for this sample is stale -- drop it rather than
    # recompute here, so it's simply recomputed, once, next time it's
    # actually needed.
    if (!is.null(rescopy$mapd) && preview$sample %in% names(rescopy$mapd))
      rescopy$mapd[[preview$sample]] <- NULL
    reslocal <<- rescopy
    trackVersion(trackVersion() + 1L)

    shinyalert("Applied",
               paste0("The refitted profile for '", preview$sample,
                      "' (binsize=", preview$params$binsize,
                      ", segmentation_alpha=", preview$params$segmentation_alpha,
                      ") has been applied to the Modifier working station (allSolutions.refitted.manual / ",
                      "allProfiles.refitted.manual). The original fit is untouched. ",
                      "Open the Modifier tab to continue refining it."),
               type="success")
  })

  #########################################################################
  ######################## SAMPLE VIEWER ##################################
  #########################################################################
  # Full-screen, one-sample-at-a-time browser with fast prev/next paging.
  # Always shows the current "manual" (working) solution/profile -- i.e.
  # whatever is presently active in the Modifier tab -- so edits made
  # there (or in Re-segment -> Apply) show up here too.
  #
  # Speed: output$viewer_profile renders to a png() device (see below for
  # why -- plotSolution()'s point overlay needs a fixed-size device) and
  # manually caches the resulting file per sample fingerprint in
  # viewerPlotCache: an index plus purity/ploidy/a numeric fingerprint of
  # the profile. Paging back to an already-visited, unedited sample is
  # then served straight from the cached PNG -- no redraw -- which is what
  # keeps this fast regardless of how many samples there are. Any edit
  # changes the fingerprint, so the cache can never show a stale profile.
  #########################################################################

  viewerIndex     <- reactiveVal(1L)
  excludedSamples <- reactiveVal(character(0))
  # Bumped whenever allTracks.processed is rewritten (a fresh load, or a
  # Re-segment/Refit "Apply") so the MAPD badge below knows to recompute
  # instead of showing a stale value for a sample whose track just changed.
  trackVersion    <- reactiveVal(0L)

  # Keeps the "Exclude this sample" checkbox showing the right state for
  # whichever sample is currently displayed. Safe to call unconditionally
  # after any navigation -- if the value doesn't actually change, no
  # extra input event fires.
  syncExcludeCheckbox <- function(samples, idx) {
    if (length(samples) == 0 || is.na(idx) || idx < 1 || idx > length(samples)) return(invisible(NULL))
    updateCheckboxInput(session, "viewer_exclude",
                         value = samples[idx] %in% isolate(excludedSamples()))
  }

  navigateViewer <- function(delta) {
    samples <- getSamples()
    n <- length(samples)
    if (n == 0) return(invisible(NULL))
    cur <- viewerIndex()
    if (is.na(cur) || cur < 1 || cur > n) cur <- 1L
    new_i <- ((cur - 1L + delta) %% n) + 1L
    viewerIndex(new_i)
    updateSelectInput(session, "viewer_sample_select", selected = samples[new_i])
    syncExcludeCheckbox(samples, new_i)
  }

  observeEvent(input$viewer_next_btn, { navigateViewer(1L) })
  observeEvent(input$viewer_prev_btn, { navigateViewer(-1L) })
  observeEvent(input$viewer_next_key, { navigateViewer(1L) })
  observeEvent(input$viewer_prev_key, { navigateViewer(-1L) })

  observeEvent(input$viewer_sample_select, {
    samples <- getSamples()
    idx <- which(samples == input$viewer_sample_select)[1]
    if (!is.na(idx)) {
      viewerIndex(idx)
      syncExcludeCheckbox(samples, idx)
    }
  }, ignoreInit = TRUE)

  # Toggling the checkbox adds/removes the *currently displayed* sample
  # from the exclusion set. This handler is idempotent by design (it
  # only adds a sample that isn't already excluded, or removes one that
  # is), so it stays correct whether it fires from a real user click or
  # from syncExcludeCheckbox() updating the widget after navigation.
  observeEvent(input$viewer_exclude, {
    samples <- getSamples()
    idx <- viewerIndex()
    if (length(samples) == 0 || is.na(idx) || idx < 1 || idx > length(samples)) return(NULL)
    samp <- samples[idx]
    cur  <- excludedSamples()
    if (isTRUE(input$viewer_exclude)) {
      if (!(samp %in% cur)) excludedSamples(c(cur, samp))
    } else {
      if (samp %in% cur) excludedSamples(setdiff(cur, samp))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$viewer_clear_exclusions, {
    excludedSamples(character(0))
    syncExcludeCheckbox(getSamples(), viewerIndex())
  })

  # `reslocal` is a plain global, not a reactiveVal, so it creates no
  # dependency on its own -- refresh the dropdown/position explicitly
  # whenever this tab is opened, same pattern as the Re-segment tab.
  observeEvent(input$nav_page, {
    if (identical(input$nav_page, "Sample Viewer")) {
      samples <- getSamples()
      n <- length(samples)
      if (n == 0) return(NULL)
      idx <- viewerIndex()
      if (is.na(idx) || idx < 1 || idx > n) { idx <- 1L; viewerIndex(1L) }
      updateSelectInput(session, "viewer_sample_select", choices = samples, selected = samples[idx])
      syncExcludeCheckbox(samples, idx)
    }
  })

  output$viewer_position_label <- renderText({
    samples <- getSamples()
    n <- length(samples)
    if (n == 0) return("No data loaded yet -- load a file from the Welcome tab.")
    idx <- viewerIndex()
    if (is.na(idx) || idx < 1 || idx > n) idx <- 1L
    paste0("Sample ", idx, " of ", n, ":  ", samples[idx])
  })

  # MAPD (Median Absolute Pairwise Difference) badge for whatever sample
  # is currently displayed -- a quick bin-to-bin noise readout to help
  # decide whether to exclude it. Colour is relative to the OTHER
  # currently loaded samples (tertiles of this dataset's own MAPD
  # distribution), not an absolute/universal cutoff, since what counts as
  # "noisy" varies by assay and by binsize/segmentation_alpha -- so this
  # flags samples that stand out within this dataset rather than
  # asserting a fixed threshold.
  output$viewer_mapd <- renderUI({
    trackVersion()  # dependency only: recompute after a track is refit
    samples <- getSamples()
    n <- length(samples)
    idx <- viewerIndex()
    if (is.null(reslocal) || n == 0 || is.na(idx) || idx < 1 || idx > n) return(NULL)

    all_mapd  <- getAllSampleMAPD()
    this_mapd <- all_mapd[[samples[idx]]]

    if (is.null(this_mapd) || is.na(this_mapd)) {
      return(tags$div("MAPD: not available for this sample",
                       style="color:#833e03; font-size:12px; font-style:italic;"))
    }

    finite_vals <- all_mapd[is.finite(all_mapd)]
    badge_color <- "#2e7d32"  # green: lower-noise third (or only sample)
    if (length(finite_vals) >= 3) {
      q <- stats::quantile(finite_vals, probs = c(1/3, 2/3), na.rm = TRUE)
      if (this_mapd > q[2])      badge_color <- "#c62828"  # red: noisier third
      else if (this_mapd > q[1]) badge_color <- "#ef6c00"  # orange: middle third
    }

    tagList(
      tags$div(style=paste0("display:inline-block; padding:6px 18px; border-radius:20px; ",
                             "background-color:", badge_color, "; color:#FFFFFF; ",
                             "font-weight:bold; font-size:14px;"),
               paste0("MAPD: ", signif(this_mapd, 3))),
      tags$div("Median absolute difference between consecutive bins on the unsegmented track -- lower is cleaner. Colour is relative to the other loaded samples.",
                style="color:#833e03; font-size:11px; max-width:520px; margin:4px auto 0 auto; line-height:1.4;")
    )
  })

  output$viewer_exclude_count <- renderText({
    n_excl <- length(excludedSamples())
    n_tot  <- length(getSamples())
    if (n_excl == 0) "No samples excluded."
    else paste0(n_excl, " of ", n_tot, " sample(s) marked for exclusion.")
  })

  # Cache of rendered PNGs for this session: cache_key -> file path. A
  # plain (non-reactive) environment, so filling it does not re-trigger
  # the render.
  viewerPlotCache <- new.env(parent = emptyenv())

  output$viewer_profile <- renderImage({
    samples <- getSamples()
    idx <- viewerIndex()
    if (is.null(reslocal) || length(samples) == 0 ||
        is.na(idx) || idx < 1 || idx > length(samples) || !isSampleIndexValid(idx)) {
      blank <- tempfile(fileext = '.png')
      png(blank, width = 10, height = 10); plot.new(); dev.off()
      return(list(src = blank, contentType = "image/png", alt = "",
                  style = "width:1px; height:1px;"))
    }

    sol <- reslocal$allSolutions.refitted.manual[[idx]]

    # Guard against reslocal$sex being shorter than the sample list (it's
    # not covered by isSampleIndexValid()) -- indexing past its end with
    # [[idx]] would otherwise throw "subscript out of bounds" here.
    ismale <- resolve_ismale(reslocal$sex, idx)

    prof <- reslocal$allProfiles.refitted.manual[[idx]]
    cache_key <- paste(idx,
                        if (!is.null(sol))  sol$purity else NA,
                        if (!is.null(sol))  sol$ploidy else NA,
                        if (!is.null(prof)) sum(suppressWarnings(as.numeric(unlist(prof))), na.rm = TRUE) else NA,
                        sep = "|")

    # The image is drawn at a fixed 950x400px and scaled to the container
    # width with CSS; the height follows to keep the aspect ratio.
    img_style <- "width:100%; height:auto; display:block; margin:0 auto;"

    cached <- get0(cache_key, envir = viewerPlotCache, inherits = FALSE)
    if (!is.null(cached) && file.exists(cached)) {
      return(list(src = cached, contentType = "image/png",
                  alt = "Sample profile", style = img_style))
    }

    # Drawn on a png() device at the same fixed 950x400 size as the
    # Modifier tab's "Original" panel: plotSolution()'s layout (including
    # the logR points) depends on the device size, so the image is scaled
    # for display in CSS (img_style) rather than drawn larger.
    outfile <- tempfile(fileext = '.png')
    png(outfile, width = 950, height = 400)
    plotSolution(reslocal$allTracks.processed[[idx]],
                 purity = sol$purity,
                 ploidy = sol$ploidy,
                 ismale = ismale,
                 gamma  = reslocal$gamma,
                 sol    = sol)
    dev.off()

    assign(cache_key, outfile, envir = viewerPlotCache)

    list(src = outfile, contentType = "image/png", alt = "Sample profile", style = img_style)
  }, deleteFile = FALSE)

  # ── Save filtered .Rda: same content the Modifier tab's "Save .Rda"
  # produces, minus whichever samples are currently checked "Exclude
  # this sample from saved output". reslocal itself is left untouched --
  # only this exported copy is filtered -- so excluded samples remain
  # fully available in the app for further editing.
  output$viewer_save <- downloadHandler(
    filename = function() { "result_manualfitting_filtered.Rda" },
    content  = function(file) {
      if (is.null(reslocal)) {
        shinyalert("Error", "Please load data first.", type = "error")
        req(FALSE)
      }
      to_exclude <- excludedSamples()
      if (length(to_exclude) > 0 && length(to_exclude) >= length(getSamples())) {
        shinyalert("Nothing to save",
                   "Every sample is currently marked for exclusion -- uncheck at least one before saving.",
                   type = "warning")
        req(FALSE)
      }
      showModal(modalDialog(div(tags$b("Loading...", style="color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      res <- excludeSamplesFromRes(reslocal, to_exclude)
      save(res, file=file)
    }
  )

}

shinyApp(ui=ui, server=server)
