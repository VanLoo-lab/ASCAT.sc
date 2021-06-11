gcTrack <-
function(chr,
                    starts,
                    ends,
                    dna,
                    window=5000)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]],
                                               window,
                                               c("G","C")))/window
    if(ends[length(ends)]>length(gc))
    {
        ends[length(ends)] <- length(gc)
    }
    gc <- sapply(1:length(starts),function(x) mean(gc[starts[x]:ends[x]]))
    names(gc) <- paste("bin",1:length(starts),sep="-")
    gc
}
