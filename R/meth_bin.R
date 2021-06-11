meth_bin <- function(logr,starts,ends,chrs,grbins)
{
    grprobes <- GenomicRanges::GRanges(chrs,IRanges(starts,ends))
    overlaps <- GenomicRanges::findOverlaps(grprobes,grbins)
    nl <- tapply(logr[queryHits(overlaps)],subjectHits(overlaps),median)
    starts <- start(grbins[as.numeric(names(nl))]@ranges)
    ends <- end(grbins[as.numeric(names(nl))]@ranges)
    chrs <- as.vector(grbins[as.numeric(names(nl))]@seqnames)
    list(nl,starts,ends,chrs)
}
