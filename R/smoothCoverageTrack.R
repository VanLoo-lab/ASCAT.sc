smoothCoverageTrack <- function(lCT,
                                lSe,
                                lGCT,
                                lNormals=NULL,
                                method=c("loess",
                                         "lowess"))
{
    allRec <- unlist(lapply(lCT,function(x) log2(x$records+1)))
    if(!is.null(lNormals))
    {
        if(!"records"%in%names(lNormals[[1]]))
            allRec <- smoothNormals(allRec, lNormals)
        else
            allRec <- allRec-unlist(lapply(lNormals,function(x) log2(x$records+1)))
    }
    allGC <- unlist(lGCT)
    starts <- c(0,cumsum(sapply(lCT,nrow)[-c(length(lCT))]))+1
    ends <- cumsum(sapply(lCT,nrow))
    if(method[1]=="lowess")
        smoothT <- lowess(allGC,allRec,f=2/3,iter=10)
    if(method[1]=="loess")
        smoothT <- myloess(LL=length(allRec),allRec~allGC)
    for(i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]],smoothT$fitted[starts[i]:ends[i]],
                          smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]])-1):ncol(lCT[[i]])] <- c("fitted","smoothed")
    }
    return(lCT)
}
