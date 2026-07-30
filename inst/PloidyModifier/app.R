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
library(spatstat.geom)
library(dipsaus)

reslocal <- NULL

pathrdata <- NULL
rdata <- NULL
coords <- NULL
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

plotScGrid <- function(solution) {
  x_values <- seq(1.1, 5, by = 0.1)
  y_values <- seq(0.95, 1, by = 0.01)
  grid <- expand.grid(x = x_values, y = y_values)
  plot(NA, xlim = c(1.1, 5), ylim = c(0.95, 1), xlab = "Ploidy",
       ylab = "Purity", type = "n")
  for (i in 1:nrow(grid)) {
    rect(grid$x[i] - 0.05, grid$y[i] - 0.005,
         grid$x[i] + 0.05, grid$y[i] + 0.005, col = "skyblue", border = "steelblue1")
  }
  purity  <- solution$purity
  ploidy  <- solution$ploidy
  points(ploidy, purity, pch = 4, col = "darkred", cex = 2, lwd = 2)
}

#######################################################################
######################### plotSunrise OVERRIDE ########################
# Defines plotSunrise locally so it shadows the ASCAT.sc package      #
# version, adding: single-purity band, small-matrix axes, null-bao    #
# guard, and withCallingHandlers so warnings don't abort the plot.     #
#######################################################################

plotSunrise <- function(solution, localMinima = FALSE, plotClust = FALSE,
                        is_sc = FALSE, N = 10)
{
  # ── Axis helpers ──────────────────────────────────────────────────────
  .ploidy_axis <- function(errs) {
    n <- ncol(errs)
    if (n <= 10L) {
      axis(side = 1, at = (seq_len(n) - 0.5) / n,
           labels = signif(as.numeric(colnames(errs)), 2))
    } else {
      idx <- pmax(1L, pmin(n, round(seq(0.1, 1, 0.1) * n)))
      axis(side = 1, at = seq(0.1, 1, 0.1),
           labels = signif(as.numeric(colnames(errs)[idx]), 2))
    }
  }
  .purity_axis_methyl <- function(errs) {
    n <- nrow(errs)
    if (n <= 10L) {
      axis(side = 2, at = 1 - (seq_len(n) - 0.5) / n,
           labels = signif(as.numeric(rownames(errs)), 2))
    } else {
      idx <- pmax(1L, pmin(n, round(seq(0.1, 1, 0.1) * n)))
      axis(side = 2, at = seq(0.1, 1, 0.1),
           labels = signif(as.numeric(rownames(errs)[(n:1)[idx]]), 2))
    }
  }
  .purity_axis_sc <- function(errs) {
    n <- nrow(errs)
    if (n < 2L) return(invisible(NULL))   # guard: seq(..., by=Inf) -> NaN
    axis(side = 2, at = seq(0.1, 1, 0.9 / (n - 1)), labels = rev(rownames(errs)))
  }

  tryCatch(
    withCallingHandlers({

      rdbu10 <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7",
                  "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")
      hmcol <- colorRampPalette(rdbu10)(256)
      hmcol[1:70]    <- colorRampPalette(hmcol[28:70])(70)
      hmcol[197:256] <- colorRampPalette(hmcol[197:232])(60)
      .getCol <- function(x) { x <- pmax(0, pmin(1, x)); hmcol[round(x * 255) + 1] }

      errs     <- solution$errs
      errs     <- errs - min(errs)
      errs.max <- max(solution$errs[!is.infinite(solution$errs)])
      errs[is.infinite(errs)] <- errs.max
      errs <- errs / errs.max

      purity_vals <- as.numeric(rownames(errs))
      if (purity_vals[1] < purity_vals[nrow(errs)])
        errs <- errs[rev(seq_len(nrow(errs))), ]

      single_purity <- nrow(errs) == 1L
      im <- matrix(.getCol(as.vector(1 - errs)), nrow(errs), ncol(errs))

      plot(0, 0, col = rgb(0, 0, 0, 0), xlab = "ploidy", ylab = "purity",
           xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, 1), ylim = c(0, 1))

      # ── Single-purity: horizontal colour band ───────────────────────
      if (single_purity) {
        BAND_LO <- 0.35; BAND_HI <- 0.65; y_mid <- 0.5
        suppressWarnings(rasterImage(as.raster(im), 0, BAND_LO, 1, BAND_HI))
        sol_col <- which(colnames(errs) == as.character(solution$ploidy))
        if (length(sol_col))
          points(sol_col / ncol(errs), y_mid, col = "chartreuse", pch = "X", cex = 1.5)
        if (localMinima) {
          row_vals <- as.numeric(errs[1, ]); n_cols <- length(row_vals)
          if (n_cols == 1L) {
            best_cols <- 1L
          } else {
            is_min <- logical(n_cols)
            is_min[1]      <- row_vals[1]      <= row_vals[2]
            is_min[n_cols] <- row_vals[n_cols] <= row_vals[n_cols - 1]
            if (n_cols > 2L)
              for (j in seq(2L, n_cols - 1L))
                is_min[j] <- row_vals[j] <= row_vals[j-1] && row_vals[j] <= row_vals[j+1]
            best_cols <- which(is_min)
            if (length(best_cols) > N)
              best_cols <- best_cols[order(row_vals[best_cols])][seq_len(N)]
          }
          bao1 <- list(bao    = cbind(rep(1L, length(best_cols)), best_cols),
                       ao     = cbind(rep(1L, length(best_cols)), best_cols),
                       clusts = rep(1L, length(best_cols)))
          text(bao1$bao[, 2] / ncol(errs), y_mid,
               labels = seq_len(nrow(bao1$bao)), col = "white")
        }
        .ploidy_axis(errs)
        axis(side = 2, at = y_mid, labels = rownames(errs)[1])

      # ── Normal multi-purity: full heatmap ───────────────────────────
      } else {
        rasterImage(as.raster(im), 0, 0, 1, 1)

        if (!is_sc) {
          sol_row <- which(rownames(errs) == as.character(solution$purity))
          sol_col <- which(colnames(errs) == as.character(solution$ploidy))
          if (length(sol_row) && length(sol_col))
            points(sol_col / ncol(errs), 1 - sol_row / nrow(errs),
                   col = "chartreuse", pch = "X", cex = 1.5)
          if (localMinima) {
            bao1 <- findLocalMinima(errs, N = N)
            if (!is.null(bao1) && !is.null(bao1$bao) && nrow(bao1$bao) > 0) {
              ao <- bao1$bao
              text(ao[, 2] / ncol(errs), 1 - ao[, 1] / nrow(errs),
                   labels = seq_len(nrow(ao)), col = "white", pch = 19)
              if (plotClust && !is.null(bao1$ao)) {
                ao <- bao1$ao
                text(ao[, 2] / ncol(errs), 1 - ao[, 1] / nrow(errs),
                     labels = seq_len(nrow(ao)),
                     col = RColorBrewer::brewer.pal(12, "Paired")[bao1$clusts], cex = 0.6)
              }
            }
          }
          .ploidy_axis(errs); .purity_axis_methyl(errs)

        } else {
          i_sol   <- which(rownames(errs) == as.character(solution$purity))
          sol_col <- which(colnames(errs) == as.character(solution$ploidy))
          if (length(i_sol) && length(sol_col)) {
            y_sol <- 1 - (i_sol - 1) * 0.9 / (nrow(errs) - 1)
            points(sol_col / ncol(errs), y_sol, col = "chartreuse", pch = "X", cex = 1.5)
          }
          if (localMinima) {
            bao1 <- findLocalMinima(errs, N = N)
            if (!is.null(bao1) && !is.null(bao1$bao) && nrow(bao1$bao) > 0) {
              ao <- bao1$bao
              text(ao[, 2] / ncol(errs),
                   1 - (ao[, 1] - 1) * 0.9 / (nrow(errs) - 1),
                   labels = seq_len(nrow(ao)), col = "white", pch = 19)
              if (plotClust && !is.null(bao1$ao)) {
                ao <- bao1$ao
                text(ao[, 2] / ncol(errs),
                     1 - (ao[, 1] - 1) * 0.9 / (nrow(errs) - 1),
                     labels = seq_len(nrow(ao)),
                     col = RColorBrewer::brewer.pal(12, "Paired")[bao1$clusts], cex = 0.6)
              }
            }
          }
          .ploidy_axis(errs); .purity_axis_sc(errs)
        }
      }

      text(0.8, 0.1, cex = 0.9,
           as.expression(bquote(paste("max ", phi[T], " hit"))), col = rgb(1, 1, 1))
      if (localMinima) return(bao1)

    }, warning = function(w) {
      message("plotSunrise warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) { print(e) }
  )
}

#######################################################################
#################### findLocalMinima OVERRIDE ########################
# Robust for small (e.g. 2x2) matrices: ensures ao is always a       #
# matrix, guards hclust against single-point input, uses drop=FALSE  #
# throughout so subsetting never silently collapses to a vector.     #
#######################################################################

findLocalMinima <- function(mat, N = 5)
{
  tryCatch({

    shifts <- list(c(1,1), c(0,1), c(1,0), c(-1,-1),
                   c(-1,+1), c(-1,0), c(0,-1), c(1,-1))
    nr <- nrow(mat)
    nc <- ncol(mat)

    applyShifts <- function(mat, shifts, nr, nc) {
      if (shifts[1] == -1) mat <- rbind(mat[-1L, ],    Inf)
      if (shifts[1] ==  1) mat <- rbind(Inf,            mat[-nr, ])
      if (shifts[2] == -1) mat <- cbind(mat[, -1L],    Inf)
      if (shifts[2] ==  1) mat <- cbind(Inf,            mat[, -nc])
      mat
    }

    nmat <- vector("list", 8L)
    for (i in 1:8) {
      nmat[[i]] <- apply(applyShifts(mat, shifts[[i]], nr, nc) > mat, 2, as.numeric)
    }

    isLocalOptima <- Reduce("*", nmat)
    ao <- which(isLocalOptima == 1, arr.ind = TRUE)

    # which() returns a *named vector*, not a matrix, when only 1 cell qualifies.
    # All subsequent [, 1] / [, 2] indexing fails on a vector.
    if (!is.matrix(ao))
      ao <- matrix(ao, nrow = 1L, dimnames = list(NULL, c("row", "col")))

    ord1  <- order(mat[isLocalOptima == 1], decreasing = FALSE)
    ao1   <- ao <- ao[ord1, , drop = FALSE]
    errs1 <- mat[ao]

    aopp <- cbind(as.numeric(rownames(mat)[ao[, 1]]) * 5,
                  as.numeric(colnames(mat)[ao[, 2]]))

    # hclust requires >= 2 points; handle single-optimum case explicitly
    if (nrow(ao) < 2L) {
      clusts1 <- 1L
      bests1  <- 1L
    } else {
      clusts1 <- cutree(hclust(dist(aopp), method = "ward.D2"), h = 0.1)
      bests1  <- unlist(
        tapply(seq_along(clusts1), clusts1,
               function(x) x[which.min(errs1[x])]),
        use.names = FALSE)
    }

    ao   <- ao[bests1, , drop = FALSE]
    ord  <- order(mat[isLocalOptima == 1][ord1][bests1], decreasing = FALSE)
    ao   <- ao[ord,   , drop = FALSE]

    aopp <- cbind(as.numeric(rownames(mat)[ao[, 1]]) * 5,
                  as.numeric(colnames(mat)[ao[, 2]]))

    if (nrow(ao) < 2L) {
      clusts <- 1L
      errs   <- mat[ao]
      bests  <- 1L
    } else {
      clusts <- cutree(hclust(dist(aopp), method = "ward.D2"), h = 0.15)
      errs   <- mat[ao]
      bests  <- unlist(
        tapply(seq_along(clusts), clusts,
               function(x) x[which.min(errs[x])]),
        use.names = FALSE)
      bests  <- bests[order(mat[ao[bests, , drop = FALSE]], decreasing = FALSE)]
    }

    if (length(bests) < N)
      bests <- c(rep(bests[1], N - length(bests)), bests)

    ao <- ao[bests, , drop = FALSE]

    list(bao    = ao[1:N, , drop = FALSE],
         ao     = ao1,
         clusts = clusts1[ord1])

  },
  error   = function(e) { print(e) },
  warning = function(w) { print(w) })
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
                              "em {display:block; padding:9.5px; margin: auto; width: 800px; height: 95px; font-size:14px; color: #833e03; line-height:5px; word-break:break-all; word-wrap:break-word; white-space:pre-wrap; background-color:#fae6d4; border:4px solid #833e03; border-radius:4px;}",".sw-dropdown-content {  max-width: min(1200px, calc(100vw - 24px)) !important;  overflow-x: hidden;}"
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
                                                  style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 2px")
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
                                       selectInput("Chr1", "Choose first chromosome", c(1:22,"X","Y"),
                                                   selected=NULL, multiple=FALSE, selectize=FALSE, width="75%"),
                                       textInput("cn1", "Choose first copy number", value="", placeholder=NULL, width="75%"),
                                       selectInput("Chr2", "Choose second chromosome", c(1:22,"X","Y"),
                                                   selected=NULL, multiple=FALSE, selectize=FALSE, width="75%"),
                                       textInput("cn2", "Choose second copy number", value="", placeholder=NULL, width="75%"),
                                       actionButton("modify", label="Apply",
                                                    style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Refit segments", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       tags$h3("Modify segment on graph"),
                                       h4("To modify the copy number of 2 different segments directly on the Working station profile, please click on the segment you wish to modify, then on its desired position. Repeat for the second segment, then click on 'Apply'"),
                                       actionButton("refit", label="Apply",
                                                    style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Refit segments on graph", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       sliderInput("ploidy", "Shift ploidy by:", -3, 4, 1, step=1,
                                                   round=FALSE, ticks=TRUE, animate=FALSE, width=NULL,
                                                   sep=",", pre=NULL, post=NULL),
                                       actionButton("shift", label="Apply",
                                                    style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                       size="lg", circle=FALSE, status="info", right=TRUE,
                                       label="Shift sample ploidy", width="450px"
                                     ),
                                     br(),
                                     dropdownButton(
                                       tags$h3("Shift ploidy on graph"),
                                       h4("To shift the ploidy of the whole sample directly on the Working station profile, please click on a point with the y axis position corresponding to the desired ploidy value, then click on 'Apply'"),
                                       actionButton("shift_graph", label="Apply",
                                                    style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
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
                                                style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px")),
                            column(width=3, offset=1,
                                   downloadButton("savetxt", label="Save profiles",
                                                  style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px")),
                            column(width=3, offset=1,
                                   downloadButton("save", label="Save .Rda",
                                                  style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"))
                          ), br(), br()
                 )
)

#######################################################################
######################### SERVER ######################################
#######################################################################

server <- function(input, output, session) {

  vals    <- reactiveVal()
  volumes <- getVolumes()
  coords  <- reactiveValues(x=NULL, y=NULL)
  output$profile  <- renderPlot(NULL)
  output$profile2 <- renderPlot(NULL)

  observeEvent(input$profile2_click, {
    coords$x <- c(coords$x, input$profile2_click$x)
    coords$y <- c(coords$y, input$profile2_click$y)
  })

  observeEvent(input$sunrise2_click, {
    coords$x <- input$sunrise2_click$x
    coords$y <- input$sunrise2_click$y
  })

  sampleName <- reactive({ input$samples })
  optValue   <- reactive({ input$optima })
  result     <- reactive({ reslocal })
  shiftv     <- reactive({ input$ploidy })
  chrs       <- reactive({ list(input$Chr1, input$Chr2, input$cn1, input$cn2) })

  # ── Data type chosen on the Welcome page ────────────────────────────
  # Values: "sc" | "methyl" | "shallow"
  # dataMode: prefers the mode embedded in the loaded object (reslocal$mode)
  # so data with an explicit mode field always routes correctly, regardless
  # of what was selected in the welcome-page dropdown.
  # Falls back to the dropdown for objects without a mode field (shallow).
  # Values: "sc" | "methyl" | "shallow"
  dataMode <- reactive({
    if (!is.null(reslocal) && "mode" %in% names(reslocal))
      reslocal$mode
    else
      input$data_type
  })

  # ── Inline Euclidean distance (replaces spatstat.geom::crossdist) ────
  # crossdist(x1,y1,x2,y2) called with scalar args is just sqrt distance;
  # replacing it avoids a "promise already under evaluation" trap that
  # spatstat.geom's generic dispatch can trigger inside tryCatch/isolate.
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
    updateNavbarPage(session=session, inputId="nav_page", selected="Modifier")
    tryCatch({

      filepath <<- as.character(pathrdata$datapath)
      load(filepath)
      reslocal <<- res

      updateSelectInput(session, "samples", label=NULL, choices=getSamples())

      # Sync the data-type dropdown with the mode stored in the file
      detected_mode <- if ("mode" %in% names(reslocal)) reslocal$mode else "shallow"
      updateSelectInput(session, "data_type", selected = detected_mode)

      mode <- dataMode()   # "sc" | "methyl" | "shallow"

      # ── Set gamma ──────────────────────────────────────────────────
      reslocal$gamma <<- if (mode == "methyl") 0.55 else 1

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

      # ── Original profile (top panel) ────────────────────────────────
      output$profile <- renderImage({
        outfile <- tempfile(fileext='.png')
        png(outfile, width=950, height=400)

        if (mode == "sc") {
          plotSolution(reslocal$allTracks.processed[[1]],
                       purity  = reslocal$allSolutions[[1]]$purity,
                       ploidy  = reslocal$allSolutions[[1]]$ploidy,
                       ismale  = if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
                       gamma   = reslocal$gamma,
                       sol     = reslocal$allSolutions[[1]])
        } else {
          # methyl / shallow: use .refitted.auto if present
          if ("allSolutions.refitted.auto" %in% names(reslocal)) {
            plotSolution(reslocal$allTracks.processed[[1]],
                         purity = reslocal$allSolutions.refitted.auto[[1]]$purity,
                         ploidy = reslocal$allSolutions.refitted.auto[[1]]$ploidy,
                         ismale = if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
                         gamma  = reslocal$gamma,
                         sol    = reslocal$allSolutions.refitted.auto[[1]])
          } else {
            plotSolution(reslocal$allTracks.processed[[1]],
                         purity = reslocal$allSolutions[[1]]$purity,
                         ploidy = reslocal$allSolutions[[1]]$ploidy,
                         ismale = if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
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
                             ismale        = if(!is.null(reslocal$sex)) reslocal$sex[[1]]=="male" else "female",
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

      removeModal()

    }, error=function(e) { })
  })

  #########################################################################
  ######################## DISCARD/KEEP ###################################
  #########################################################################

  observeEvent(input$discard, {
    vals(input$samples)
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      mode  <- dataMode()

      if (mode == "sc") {
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity        = reslocal$allSolutions[[index]]$purity,
                               ploidy        = reslocal$allSolutions[[index]]$ploidy,
                               ismale        = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                               gamma         = reslocal$gamma,
                               sol           = reslocal$allSolutions[[index]],
                               ambiguousFlag = FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({ localOpt <<- isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE)) })
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
                               ismale        = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
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
    vals(input$samples)
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      cnHidden(input$hideCN)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=input$hideCN, transparentCN=cnTrans(), zoomPoints=pSize()))
      })
    } else { return(NULL) }
  })

  observeEvent(input$transparentCN, {
    vals(input$samples)
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      cnTrans(input$transparentCN)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=cnHidden(), transparentCN=input$transparentCN, zoomPoints=pSize()))
      })
    } else { return(NULL) }
  })

  observeEvent(input$sizeP, {
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      pSize(input$sizeP)
      output$profile2 <- renderPlot({
        isolate(plotSolution(tracksSingle   = reslocal$allTracks.processed[[index]],
                             purity         = reslocal$allSolutions.refitted.manual[[index]]$purity,
                             ploidy         = reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                             ismale         = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
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
    mode  <- dataMode()

    tryCatch({
      solution <- reslocal$allSolutions[[index]]
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
                     purity, ploidy, ismale=reslocal$sex[index]=="male", gamma=reslocal$gamma),
          CHRS=reslocal$chr)
        reslocal$allSolutions.refitted.manual[[index]]$ploidy <<- ploidy
        reslocal$allSolutions.refitted.manual[[index]]$purity <<- purity
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=purity, ploidy=ploidy, gamma=reslocal$gamma,
                               ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                               sol=reslocal$allSolutions.refitted.manual[[index]],
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], is_sc=TRUE, localMinima=TRUE))
        })

      } else {
        # methyl / shallow share the same coordinate mapping
        errs <- errs[rev(seq_len(nrow(errs))), ]

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
                     purity, ploidy, ismale=reslocal$sex[index]=="male", gamma=reslocal$gamma),
          CHRS=reslocal$chr)
        reslocal$allSolutions.refitted.manual[[index]]$ploidy <<- ploidy
        reslocal$allSolutions.refitted.manual[[index]]$purity <<- purity
        output$profile2 <- renderPlot({
          isolate(plotSolution(reslocal$allTracks.processed[[index]],
                               purity=purity, ploidy=ploidy, gamma=reslocal$gamma,
                               ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
        })
        output$sunrise2 <- renderPlot({
          isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
        })
      }

      coords$y <<- NULL
      coords$x <<- NULL

    },
    error=function(e) {
      message('An Error Occurred'); print(e)
      shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type="error")
    },
    warning=function(w) {
      message('A Warning Occurred'); print(w)
      shinyalert("Warning", "New solution is ambiguous: reverted to old one", type="error")
    })
  })

  #########################################################################
  ######################## MODIFY SEGMENTS ON GRAPH #######################
  #########################################################################

  observeEvent(input$refit, {
    vals <- as.numeric(chrs())
    if (!is.null(chrs())) {
      index  <- getIndex(sampleName())
      mode   <- dataMode()
      breaks <- c(0, cumsum(sapply(reslocal$allTracks.processed[[1]]$lSegs,
                                   function(x) max(x$output$loc.end)) / 1e+06))

      tryCatch({
        chr1 <- NULL; chr2 <- NULL
        for (i in 1:length(breaks)) { if (coords$x[1] < breaks[i]) { chr1 <- i-1; break } }
        for (i in 1:length(breaks)) { if (coords$x[3] < breaks[i]) { chr2 <- i-1; break } }
        y2 <- round(coords$y[2], digits=0)
        y4 <- round(coords$y[4], digits=0)

        if (mode == "sc") {
          rescopy <- reslocal
          rescopy$allSolutions <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                          chr1=chr1, ind1=NA, n1=y2,
                                          chr2=chr2, ind2=NA, n2=y4,
                                          CHRS=c(1:22,"X","Y"), outdir="./www",
                                          gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions <- reslocal$allSolutions
          rescopy$allProfiles  <- reslocal$allProfiles
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
          })

        } else {
          # methyl / shallow
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                          chr1=chr1, ind1=NA, n1=y2,
                                          chr2=chr2, ind2=NA, n2=y4,
                                          CHRS=c(1:22,"X","Y"), outdir="./www",
                                          gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
          })
        }

        coords$y <<- NULL; coords$x <<- NULL

      },
      error=function(e) {
        message('An Error Occurred'); print(e)
        shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type="error")
      },
      warning=function(w) {
        message('A Warning Occurred'); print(w)
        shinyalert("Warning", "New solution is ambiguous: reverted to old one", type="error")
      })
    } else { return(NULL) }
  })

  #########################################################################
  ######################## MODIFY SEGMENTS ################################
  #########################################################################

  observeEvent(input$modify, {
    vals <- as.numeric(chrs())
    if (!is.null(chrs())) {
      index <- getIndex(sampleName())
      mode  <- dataMode()

      tryCatch({
        if (mode == "sc") {
          rescopy <- reslocal
          rescopy$allSolutions <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                          chr1=vals[[1]], ind1=NA, n1=vals[[3]],
                                          chr2=vals[[2]], ind2=NA, n2=vals[[4]],
                                          CHRS=c(1:22,"X","Y"), outdir="./www",
                                          gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions <- reslocal$allSolutions
          rescopy$allProfiles  <- reslocal$allProfiles
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], is_sc=TRUE))
          })

        } else {
          # methyl / shallow
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile(rescopy, sample_indice=index,
                                          chr1=vals[[1]], ind1=NA, n1=vals[[3]],
                                          chr2=vals[[2]], ind2=NA, n2=vals[[4]],
                                          CHRS=c(1:22,"X","Y"), outdir="./www",
                                          gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
          })
        }

      },
      error=function(e) {
        message('An Error Occurred'); print(e)
        shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type="error")
      },
      warning=function(w) {
        message('A Warning Occurred'); print(w)
        shinyalert("Warning", "New solution is ambiguous: reverted to old one", type="error")
      })
    } else { return(NULL) }
  })

  #########################################################################
  ######################## SHIFT ON GRAPH #################################
  #########################################################################

  observeEvent(input$shift_graph, {
    if (!is.null(chrs())) {
      tryCatch({
        index <- getIndex(sampleName())
        mode  <- dataMode()
        ypos  <- as.numeric(round(coords$y[1], digits=0))
        ploidy_cur <- round(as.numeric(reslocal$allSolutions.refitted.manual[[index]]$ploidy), digits=0)
        shiftp <- ypos - ploidy_cur

        if (mode == "sc") {
          rescopy <- reslocal
          rescopy$allSolutions <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                                CHRS=c(1:22,"X","Y"), outdir="./www",
                                                gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions <- reslocal$allSolutions
          rescopy$allProfiles  <- reslocal$allProfiles
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], is_sc=TRUE))
          })

        } else {
          # methyl / shallow
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                                CHRS=c(1:22,"X","Y"), outdir="./www",
                                                gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            isolate(plotSolution(reslocal$allTracks.processed[[index]],
                                 purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                                 ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                                 ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                                 gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                                 ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
          })
          output$sunrise2 <- renderPlot({
            isolate(localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE))
          })
        }

        coords$y <<- NULL; coords$x <<- NULL

      },
      error=function(e) {
        message('An Error Occurred'); print(e)
        shinyalert("Error", "Cannot fit profile. Please choose different values", type="error")
      },
      warning=function(w) {
        message('A Warning Occurred'); print(w)
        shinyalert("Warning", "New solution is ambiguous: reverted to old one", type="error")
      })
    } else { return(NULL) }
  })

  #########################################################################
  ######################## SHIFT ##########################################
  #########################################################################

  observeEvent(input$shift, {
    if (!is.null(chrs())) {
      index  <- getIndex(sampleName())
      mode   <- dataMode()
      shiftp <- as.numeric(shiftv())

      withCallingHandlers({
        if (mode == "sc") {
          rescopy <- reslocal
          rescopy$allSolutions <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                                CHRS=c(1:22,"X","Y"), outdir="./www",
                                                gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions <- reslocal$allSolutions
          rescopy$allProfiles  <- reslocal$allProfiles
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                         ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                         ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                         gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                         ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize())
          })
          output$sunrise2 <- renderPlot({
            isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE, is_sc=TRUE))
          })

        } else {
          # methyl / shallow
          rescopy <- reslocal
          rescopy$allSolutions.refitted.auto <- rescopy$allSolutions.refitted.manual
          rescopy$allProfiles.refitted.auto  <- rescopy$allProfiles.refitted.manual
          rescopy <- run_any_refitProfile_shift(rescopy, sample_indice=index, shift=shiftp,
                                                CHRS=c(1:22,"X","Y"), outdir="./www",
                                                gridpur=seq(-.05,.05,.01), gridpl=seq(-.1,.2,.01))
          rescopy$allSolutions.refitted.auto <- reslocal$allSolutions.refitted.auto
          rescopy$allProfiles.refitted.auto  <- reslocal$allProfiles.refitted.auto
          reslocal <<- rescopy
          output$profile2 <- renderPlot({
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions.refitted.manual[[index]]$purity,
                         ploidy=reslocal$allSolutions.refitted.manual[[index]]$ploidy,
                         ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                         gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.manual[[index]],
                         ambiguousFlag=FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize())
          })
          output$sunrise2 <- renderPlot({
            localOpt <<- plotSunrise(reslocal$allSolutions.refitted.manual[[index]], localMinima=TRUE)
          })
        }

      },
      error=function(e) {
        shinyalert("Error", "Cannot fit profile. Please choose different values", type="error")
      },
      warning=function(w) {
        shinyalert("Warning", "Cannot fit profile. Please choose different values", type="error")
      })
    } else { return(NULL) }
  })

  #########################################################################
  ######################## VIEW ###########################################
  #########################################################################

  observeEvent(input$view, {
    vals(input$samples)
    if (!is.null(sampleName())) {
      index <- getIndex(sampleName())
      mode  <- dataMode()

      # Guard: allSolutions may be shorter than allTracks.processed if some
      # samples had no valid solution in the original run.
      sol_len <- length(reslocal$allSolutions)
      man_len <- length(reslocal$allSolutions.refitted.manual)
      if (is.na(index) || index > sol_len || index > man_len) {
        shinyalert("Warning",
                   paste0("Sample '", sampleName(), "' has no solution entry ",
                          "(index ", index, " out of ", sol_len, "). ",
                          "Try a different sample."),
                   type = "warning")
        return(NULL)
      }

      output$profile <- renderImage({
        outfile <- tempfile(fileext='.png')
        png(outfile, width=950, height=400)

        if (mode == "sc") {
          plotSolution(reslocal$allTracks.processed[[index]],
                       purity=reslocal$allSolutions[[index]]$purity,
                       ploidy=reslocal$allSolutions[[index]]$ploidy,
                       ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                       gamma=reslocal$gamma, sol=reslocal$allSolutions[[index]])
        } else {
          # methyl / shallow
          if ("allSolutions.refitted.auto" %in% names(reslocal)) {
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions.refitted.auto[[index]]$purity,
                         ploidy=reslocal$allSolutions.refitted.auto[[index]]$ploidy,
                         ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                         gamma=reslocal$gamma, sol=reslocal$allSolutions.refitted.auto[[index]])
          } else {
            plotSolution(reslocal$allTracks.processed[[index]],
                         purity=reslocal$allSolutions[[index]]$purity,
                         ploidy=reslocal$allSolutions[[index]]$ploidy,
                         ismale=if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
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
                             ismale         = if(!is.null(reslocal$sex)) reslocal$sex[[index]]=="male" else "female",
                             gamma          = reslocal$gamma,
                             sol            = reslocal$allSolutions.refitted.manual[[index]],
                             ambiguousFlag  = FALSE, hideCN=cnHidden(), transparentCN=cnTrans(), zoomPoints=pSize()))
      })

      if (mode == "sc") {
        output$sunrise1 <- renderPlot({ isolate(plotSunrise(reslocal$allSolutions[[index]], is_sc=TRUE)) })
        output$sunrise2 <- renderPlot({ localOpt <<- isolate(plotSunrise(reslocal$allSolutions.refitted.manual[[index]], is_sc=TRUE, localMinima=TRUE)) })
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

}

shinyApp(ui=ui, server=server)
