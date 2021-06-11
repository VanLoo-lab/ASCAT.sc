sc_plotGenome <-
function (mat,
          colSeg = rgb(0.5, 0.2, 0.5, 0.7),
          lwdSeg = 2,
          scaleY = 0.6,
          centromeres,
                           ...)
{
    if(!is.null(centromeres))
        centromeres[, 1] <- gsub("chr", "", centromeres[, 1])
    matPos <- mat
    matPos[matPos < 0] <- 0
    matPos[matPos > 0] <- 1
    matNeg <- mat
    matNeg[matNeg > 0] <- 0
    matNeg[matNeg < 0] <- -1
    vecPositive <- rowSums(matPos)
    names(vecPositive) <- rownames(mat)
    vecNegative <- rowSums(matNeg)
    names(vecNegative) <- rownames(mat)
    par(mar = c(4, 4, 1, 1))
    chrom <- as.numeric(gsub("(.*):(.*)-(.*)", "\\1", names(vecPositive)))
    starts <- as.numeric(gsub("(.*):(.*)-(.*)", "\\2", names(vecPositive)))
    ends <- as.numeric(gsub("(.*):(.*)-(.*)", "\\3", names(vecPositive)))
    breaks <- c(0, cumsum(sapply(unique(chrom), function(x) max(ends[chrom ==
    x]))))/1e+06
    names(breaks) <- c(0, unique(chrom))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
    xlim = c(0, max(breaks)), ylim = c(-ncol(mat), ncol(mat)) *
    scaleY, xlab = "Cumulative Genomic Position (Mb)",
    ylab = "Number of Samples", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    for (i in 1:length(unique(chrom))) {
        takeChrom <- chrom == unique(chrom)[i]
        segments(starts[takeChrom]/1e+06 + breaks[i], vecPositive[takeChrom],
        ends[takeChrom]/1e+06 + breaks[i], vecPositive[takeChrom],
        col = rgb(sc_getCols(1, "R"), sc_getCols(1, "G"), sc_getCols(1,
        "B")), lwd = lwdSeg)
        segments(starts[takeChrom]/1e+06 + breaks[i], vecNegative[takeChrom],
        ends[takeChrom]/1e+06 + breaks[i], vecNegative[takeChrom],
        col = rgb(sc_getCols(-1, "R"), sc_getCols(-1, "G"), sc_getCols(-1,
        "B")), lwd = lwdSeg)
        if(!is.null(centromeres))
            abline(v = centromeres[centromeres[, 1] == i, 2]/1e+06 +
        breaks[i], lty = 2, col = rgb(0.4, 0.4, 0.4, 0.4))
        if(!is.null(centromeres))
            abline(v = centromeres[centromeres[, 1] == i, 3]/1e+06 +
        breaks[i], lty = 2, col = rgb(0.4, 0.4, 0.4, 0.4))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 1, col = rgb(0.2,
    0.1, 0.1, 1))
    text(x = breaks[2:length(breaks)] - 25, y = ncol(mat) * scaleY,
    names(breaks)[2:length(breaks)], cex = 0.4)
    max(breaks)
}
