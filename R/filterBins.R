filterBins <- function(allTracks=NULL, logr=NULL, lSe, fraction=.1, k=11)
{
    smoothlogr <- function(vec, k)
    {
        meds <- stats::runmed(vec, k=k, endrule="median")
        mads <-stats::runmed(abs(vec-meds),k=k, endrule="median")
        (vec-meds)/mads
    }
    if(!is.null(allTracks))
    {
        rr <- do.call("cbind",lapply(allTracks, function(x) unlist(lapply(x$lCTS,function(y) y$records))))
    }
    if(!is.null(logr))
    {
        rr <- logr
    }
    rr <- smoothlogr(rr, k=k)
    rranks <- apply(rr,2,rank)
    rranks.n <- rranks
    NN <- nrow(rranks)
    for(i in 1:ncol(rranks)) rranks.n[,i] <- rranks.n[,i]/(NN+1-rranks.n[,i])
    rs <- rowSums(log(rranks.n))
    llS <- sapply(1:length(lSe),function(x) length(lSe[[x]]$starts))
    lInds <- lapply(1:length(lSe),function(x) if(x==1) 1:length(lSe[[x]]$starts) else 1:length(lSe[[x]]$starts)+sum(llS[1:(x-1)]))
    threshold_low <- quantile(rs,prob=fraction/2)
    threshold_high <- quantile(rs,prob=1-fraction/2)
    lInds <- lapply(1:length(lSe),function(x)
    {
        cond <- rs[lInds[[x]]]< threshold_low | rs[lInds[[x]]]>threshold_high
        which(cond)
    })
    lInds
}
