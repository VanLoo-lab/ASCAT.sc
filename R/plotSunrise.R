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
