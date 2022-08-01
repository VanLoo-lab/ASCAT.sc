gcTrack <- function (chr,
                     starts,
                     ends,
                     dna,
                     window = 50,
                     starts.exclude=NULL,
                     ends.exclude=NULL)
{
    if(is.null(starts.exclude))
    {
        gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]], window,
                                                   c("G", "C")))/window
        if (ends[length(ends)] > length(gc)) {
            ends[length(ends)] <- length(gc)
        }
        gc <- sapply(1:length(starts), function(x) mean(gc[starts[x]:ends[x]]))
        names(gc) <- paste("bin", 1:length(starts), sep = "-")
        return(gc)
    }
    else
    {
        gc_exclude <- gcTrack(chr,starts.exclude, ends.exclude, dna, window=window)
        widths_exclude <- ends.exclude-starts.exclude
        gc <- gcTrack(chr,starts, ends, dna, window=window)
        widths <- ends-starts
        grE <- GRanges(chr, IRanges(starts.exclude, ends.exclude))
        gr <- GRanges(chr, IRanges(starts, ends))
        ovs <- findOverlaps(gr,grE)
        sH <- subjectHits(ovs)
        qH <- queryHits(ovs)
        for(i in unique(qH))
        {
            ind <- qH==i
            width_remaining <- widths[i]-sum(widths_exclude[sH[ind]])
            gc[i] <- gc[i]*widths[i]-sum(gc_exclude[sH[ind]]*widths_exclude[sH[ind]])
            gc[i] <- gc[i]/width_remaining
        }
        return(gc)
    }
}
