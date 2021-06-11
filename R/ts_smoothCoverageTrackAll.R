ts_smoothCoverageTrackAll <- function (lCT, lSe, lGCT, lNormal, lRepli=NULL)
{
    allRec <- unlist(lapply(lCT, function(x) log2(pmax(x$records,0) + 1)))
    allGC <- unlist(lGCT)
    allRepli <- unlist(lRepli)
    if(!is.null(lNormal))
        lNormal <- unlist(lapply(lNormal, function(x) log2(pmax(x$records,0) + 1)))
    starts <- c(0, cumsum(sapply(lCT, nrow)[-c(length(lCT))])) + 1
    ends <- cumsum(sapply(lCT, nrow))
    smoothT <- ts_smoothTrack(R=allRec,GC=allGC,REPLI=allRepli, NORMAL=lNormal)
    for (i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]],
                          smoothT$fitted[starts[i]:ends[i]],
                          smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]]) - 1):ncol(lCT[[i]])] <- c("fitted",
                                                                     "smoothed")
    }
    return(lCT)
}
