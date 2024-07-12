segmentTrack_pcf <- function (covtrack,
                              chr,
                              starts,
                              ends = NA,
                              sd = 0,
                              transform=FALSE,
                              ismedian=FALSE,
                              min.width = 5,
                              ALPHA=0.01,
                              SBDRY=NULL,
                              smooth=TRUE)
{
    suppressPackageStartupMessages(require(copynumber))
    covtrack <- covtrack * rnorm(length(covtrack),
                                 mean = 1,
                                 sd = sd)
    if(transform)
        covtrack <- sqrt(2^covtrack +3/8)
    x <- pcf(data.frame(chr=rep(chr,length(starts)),
                        positions=starts,
                        sample1=covtrack),
             gamma=1/ALPHA,
             verbose=F,
             kmin=min.width)
    cna <- list(data=NULL,
                output=data.frame(ID="Sample.ID",
                                  chrom=x[,"chrom"],
                                  loc.start=as.numeric(as.character(x[,"start.pos"])),
                                  loc.end=as.numeric(as.character(x[,"end.pos"])),
                                  num.mark=as.numeric(as.character(x[,"n.probes"])),
                                  seg.mean=as.numeric(as.character(x[,"mean"]))))
    if(ismedian)
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
                offset <- sum(cna.[1:(x-1),"num.mark"])
                return((offset+1):(offset+cna.[x,"num.mark"]))
            }
        })
        for(i in 1:nrow(cna.))
        {
            cna$output[i,"seg.mean"] <- median(covtrack[indices[[i]]],na.rm=T)
        }
    }
    cna
}
