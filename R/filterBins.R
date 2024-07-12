filterBins <- function(allTracks=NULL, logr=NULL, lSe, IQRC=1.5)
{
    smoothlogr <- function(vec)
    {
        meds <- stats::runmed(vec, k=5, endrule="median")
        mads <-stats::runmed(abs(vec-meds),k=5, endrule="median")
        (vec-meds)/mads
    }
    identify_IQRC <- function(data, target_fraction = 0.1)
    {
        Q1 <- quantile(data, 0.25)
        Q3 <- quantile(data, 0.75)
        IQR <- Q3 - Q1
        fraction_outside <- function(IQRC)
        {
            lower_bound <- Q1 - IQRC * IQR
            upper_bound <- Q3 + IQRC * IQR
            fraction <- mean(data < lower_bound | data > upper_bound)
            return(fraction)
        }
        search_IQRC <- function()
        {
            best_IQRC <- NULL
            best_difference <- Inf
            for (IQRC in seq(1.5, 5, by = 0.01))
            {
                fraction <- fraction_outside(IQRC)
                difference <- abs(fraction - target_fraction)
                if (difference < best_difference)
                {
                    best_difference <- difference
                    best_IQRC <- IQRC
                }
            }
            return(best_IQRC)
        }
        optimal_IQRC <- search_IQRC()
        return(optimal_IQRC)
    }
    if(!is.null(allTracks))
    {
        rr <- do.call("cbind",lapply(allTracks, function(x) unlist(lapply(x$lCTS,function(y) y$records))))
    }
    if(!is.null(logr))
    {
        rr <- logr
    }
    rr <- smoothlogr(rr)
    rranks <- apply(rr,2,rank)
    rranks.n <- rranks
    NN <- nrow(rranks)
    for(i in 1:ncol(rranks)) rranks.n[,i] <- rranks.n[,i]/(NN+1-rranks.n[,i])
    rs <- rowSums(log(rranks.n))
    llS <- sapply(1:length(lSe),function(x) length(lSe[[x]]$starts))
    lInds <- lapply(1:length(lSe),function(x) if(x==1) 1:length(lSe[[x]]$starts) else 1:length(lSe[[x]]$starts)+sum(llS[1:(x-1)]))
    IQR <- quantile(rs,prob=c(.25,.75))
    threshold_low <- IQR[1]-IQRC*diff(IQR)
    threshold_high <- IQR[2]+IQRC*diff(IQR)
    if(mean(rs<threshold_low | rs>threshold_high)>.1)
    {
        IQRC <- identify_IQRC(rs, target_fraction = .1)
        threshold_low <- IQR[1]-IQRC*diff(IQR)
        threshold_high <- IQR[2]+IQRC*diff(IQR)
    }
    if(F)
        lSe <- lapply(1:length(lSe),function(x)
        {
            cond <- rs[lInds[[x]]]< threshold_low | rs[lInds[[x]]]>threshold_high
            list(starts=lSe[[x]]$starts[!cond],
                 ends=lSe[[x]]$ends[!cond])
        })
    lInds <- lapply(1:length(lSe),function(x)
    {
        cond <- rs[lInds[[x]]]< threshold_low | rs[lInds[[x]]]>threshold_high
        which(cond)
    })
    lInds
}
