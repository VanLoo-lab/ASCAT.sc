treatGCT <- function(lGCT, window=ceiling(as.numeric(WINDOW)/as.numeric("10000")))
{
    nlGCT <- lapply(lGCT,function(gc)
    {
        l <- length(gc)
        starts <- getstartends(end=l,window=window)$starts
        ends <- getstartends(end=l,window=window)$ends
        gc <- sapply(1:length(starts),function(x) mean(gc[starts[x]:ends[x]]))
    })
    names(nlGCT) <- names(lGCT)
    nlGCT
}
