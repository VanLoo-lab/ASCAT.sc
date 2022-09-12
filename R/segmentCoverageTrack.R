segmentCoverageTrack <-
function(cT,
                                 starts,
                                 ends,
                                 chr)
{
    suppressPackageStartupMessages(require(DNAcopy))
    segments<- segmentTrack(counts,
                            chr,
                            starts,
                            ends)
    segs <- segments$output
    return(segs)
}
