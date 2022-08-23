getFilters <- function(res,
                       probs=.1,
                       outdir="./",
                       projectname="")
{
    filters <- NULL
    try({
        allT <- res$allTracks.processed
        allS <- res$allSolutions
        ## ##################################################################
        getloess <- function(qu, nr)
        {
            nms <- paste0("n",1:length(qu))
            names(qu) <- nms
            quo <- qu[order(nr,decreasing=F)]
            fitted <- stats::runmed(quo, k=31, endrule="keep")
            names(fitted) <- names(quo)
            list(fitted=fitted[nms],
                 residuals=quo[nms]-fitted[nms])
        }
        getQuality.SD <- function(allT)
        {
            sapply(allT,function(x)
            {
                median(abs(diff(unlist(lapply(x$lCTS,function(y) y$smoothed)))))
            })
        }
        ## ##################################################################
        nrecords  <-  sapply(allT,function(x) sum(unlist(lapply(x$lCTS,function(y) y$records))))
        thresholdNrec <- quantile(nrecords,probs=probs)## removing 10% cells with lowest read counts
        ambiguous <- sapply(allS,function(x) x$ambiguous)
        doublet <- sapply(allS,function(x) if(!is.null(x$bestfit)) !x$bestfit$ambiguous else F)
        qualities <- getQuality.SD(allT)
        thresholdQual <- quantile(qualities,probs=1-probs)## removing 10% cells with lowest read counts
        keep <- qualities<=thresholdQual & nrecords>=thresholdNrec & !ambiguous & !doublet
        keep2 <- !(qualities<thresholdQual & nrecords<thresholdNrec)
        keep2 <- keep2 & !ambiguous
        ll <- getloess(qualities[keep2], log2(nrecords)[keep2])
        ord <- order(log2(nrecords)[keep2],decreasing=F)
        filters <- (1:length(nrecords))%in%(which(keep2)[abs(ll$residuals)<=0.02]) & keep
        ## ##################################################################
        pdf(paste0(outdir,"/",projectname,"_sc_filters.pdf"))
        plot(qualities,
             log2(nrecords),
             xlab="Noise logr",
             ylab="Total number of reads",
             pch=ifelse(doublet,15,19),
             cex=ifelse(doublet,.3,.1),
             col=ifelse(ambiguous,rgb(1,0,.5,.5),rgb(0,0,0,.5)))
        points(ll$fitted[ord],log2(nrecords)[keep2][ord],type="l",col=rgb(1,0,0,.5),lwd=1.5)
        abline(v=thresholdQual,h=log2(thresholdNrec))
        points(qualities[keep2],
               log2(nrecords)[keep2],
               col=ifelse(abs(ll$residuals)>0.02, rgb(1,0,0,1), rgb(0,0,0,0)),cex=.15)
        try({plot(density(log2(nrecords[!doublet])))
            polygon(density(log2(nrecords[doublet])),border="red")
        },silent=T)
        plot(qualities,
             log2(nrecords),
             xlab="Noise logr",
             ylab="Total number of reads",
             pch=ifelse(filters,15,19),
             cex=ifelse(filters,.3,.1),
             col=ifelse(filters,rgb(1,0,.5,.5),rgb(0,0,0,.5)))
        dev.off()
        ## ##################################################################
    })
    res$filters <- filters
    res
}
