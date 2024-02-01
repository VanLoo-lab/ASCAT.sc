plotSolution <- function(tracksSingle,
                         purity,
                         ploidy,
                         sol=NULL,
                         lwdSeg=2,
                         ylim=c(0,8),
                         gamma=1,
                         ismale=F,
                         isPON=F,
                         allchr=NULL,
                         colFit=rgb(0.756, 0.494, 0.756),
                         col2=rgb(0.2, 0.2, 0.2, 0.2),
                         colLine= rgb(0.6, 0.6, 0.6, 0.4),
                         ...)
{
    meansSeg <- fitProfile(tracksSingle,purity,ploidy,gamma=gamma, ismale=ismale, isPON=isPON)
    tracksSingle <- normaliseByPloidy(tracksSingle)
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
         xlim = c(0, max(breaks)), xlab = "Genomic Position",
         ylab = "Total copy number", frame = F, ylim=ylim, ...)
    ##axis(side = 1)
    axis(side = 2, las=2)
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
                 col = col2)
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
           col = colLine)
    abline(h=1:50, lwd = 1,
           lty = 2,
           col = colLine)
    labels <- allchr
    if(is.null(allchr))
        labels <- if(is.null(names(tracksSingle$lCTS))) names(breaks)[2:length(breaks)]
                  else names(tracksSingle$lCTS)
    text(x = breaks[2:length(breaks)] - 25,
         y = max(ylim),
         ##if(is.null(allchr)) names(breaks)[2:length(breaks)] else allchr,
         labels=labels,
         cex = 0.6)
    dpb <- median(unlist(lapply(tracksSingle$lCTS,function(x) x$records)),na.rm=T)
    dpb <- if(all(tracksSingle$lCTS[[1]]$records==tracksSingle$lCTS[[1]]$smoothed)) NA else dpb
    mtext(side = 3,
          paste0("purity=",
                 signif(purity,2),
                 "; average ploidy=",
                 signif(ploidy,2),
                 "; tumor ploidy=",
                 signif(getTumourPhi(ploidy,purity),2),
                 if(!is.null(sol)) paste0("; ",ifelse(is.null(sol$ambiguous),"",paste0("ambiguous:",sol$ambiguous))),
                 if(is.na(dpb)) "" else paste0("; dpb=",signif(dpb,2))),
          cex=.8)
}
