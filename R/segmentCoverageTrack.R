segmentCoverageTrack <-
function(cT,
                                 starts,
                                 ends,
                                 chr)
{
    require(DNAcopy)
    segments<- segmentTrack(counts,
                            chr,
                            starts,
                            ends)
    segs <- segments$output
    return(segs)
}
