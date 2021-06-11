plotSolution <- function(tracksSingle,
                         purity,
                         ploidy,
                         lwdSeg=2,
                         ylim=c(0,8),
                         gamma=1,
                         ismale=F,
                         isPON=F,
                         colFit=rgb(0.756, 0.494, 0.756),
                         ...)
{
    meansSeg <- fitProfile(tracksSingle,purity,ploidy,gamma=gamma, ismale=ismale, isPON=isPON)
    tracksSingle <- normaliseByPloidy(tracksSingle)
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
         xlim = c(0, max(breaks)), xlab = "Genomic Position",
         ylab = "Total copy number", frame = F, ylim=ylim, ...)
    ##axis(side = 1)
    axis(side = 2)
    for (i in 1:length(tracksSingle$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
                 transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                                       purity,
                                       ploidy,
                                       gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23),
                 tracksSingle$lCTS[[i]]$end/1e+06 + breaks[i],
                 transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                                       purity,
                                       ploidy,
                                       gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23),
                 col = rgb(0.7, 0.7, 0.7, 0.6))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 + breaks[i],
                 transform_bulk2tumour(sapply(meansSeg[[i]], function(x) x$mu),
                                       purity,
                                       ploidy,
                                       gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23),
                 tracksSingle$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
                 transform_bulk2tumour(sapply(meansSeg[[i]], function(x) x$mu),
                                       purity,
                                       ploidy,
                                       gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23),
                 lwd = lwdSeg,
                 col = rgb(0.4, 0.4, 0.4, 0.4))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 + breaks[i],
                 round(sapply(meansSeg[[i]], function(x) x$roundmu)),
                 tracksSingle$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
                 round(sapply(meansSeg[[i]], function(x) x$roundmu)),
                 lwd = 2.5,
                 col = colFit)
    }
    abline(v = breaks,
           lwd = 1,
           lty = 2,
           col = rgb(0.6, 0.6, 0.6, 0.4))
    abline(h=1:50, lwd = 1,
           lty = 2,
           col = rgb(0.6, 0.6, 0.6, 0.2))
    text(x = breaks[2:length(breaks)] - 25,
         y = max(ylim),
         names(breaks)[2:length(breaks)],
         cex = 0.4)
    mtext(side = 3,
          paste0("purity=",
                 signif(purity,2),
                 "; average ploidy=",
                 signif(ploidy,2),
                 "; tumor ploidy=",
                 signif(getTumourPhi(ploidy,purity),2)))
}
