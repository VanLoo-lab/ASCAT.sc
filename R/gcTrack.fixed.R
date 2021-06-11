gcTrack.fixed <- function(chr,
                          dna,
                          step=500000)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]],
                                               step,
                                               c("G","C")))/step
    lGC <- length(gc)
    gc <- gc[seq(step,lGC,step)-round(step/2)]
    names(gc) <- paste("bin",seq(step,lGC,step)-step+1,sep="-")
    gc
}
