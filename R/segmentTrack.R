segmentTrack <- function (covtrack,
                          chr,
                          starts,
                          ends = NA,
                          sd = 0,
                          transform=FALSE,
                          min.width = 5,
                          ALPHA=0.01,
                          SBDRY=NULL,
                          smooth=TRUE)
{
    require(DNAcopy)
    covtrack <- covtrack * rnorm(length(covtrack),
                                 mean = 1,
                                 sd = sd)
    if(!transform)
        cna <- CNA(covtrack,
                   chr = rep(chr, length(covtrack)),
                   maploc = starts,
                   data.type = "logratio")
    if(transform)
        cna <- CNA(sqrt(2^covtrack +3/8),
                   chr = rep(chr, length(covtrack)),
                   maploc = starts,
                   data.type = "logratio")
    if(smooth)
        cna <- smooth.CNA(cna)
    if(is.null(SBDRY))
        capture.output(cna <- segment(cna, min.width = min.width,alpha=ALPHA))
    else
        capture.output(cna <- segment(cna, min.width = min.width,alpha=ALPHA, sbdry=SBDRY))
    if(transform)
    {
        cna. <- cna$output
        indices <- lapply(1:nrow(cna.),function(x)
        {
            if(x==1)
            {
                return(1:cna.[x,"num.mark"])
            }
            else
            {
                offset <- cumsum(cna.[1:(x-1),"num.mark"])
                return((offset+1):(offset+cna.[x,"num.mark"]))
            }
        })
        for(i in 1:nrow(cna.))
        {
            cna$output[i,"seg.mean"] <- mean(covtrack[indices[[i]]],na.rm=T)
        }
    }
    cna
}
