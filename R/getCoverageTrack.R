getCoverageTrack <-
function (bamPath,
                              chr,
                              starts,
                              ends,
                              CHRSTRING = "",
                              isDuplicate=F,
                              isSecondaryAlignment=F,
                              isNotPassingQualityControls=NA,
                              isUnmappedQuery=NA,
                              mapqFilter=30)
{
    require(Rsamtools)
    require(GenomicRanges)
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate=isDuplicate,
                                           isSecondaryAlignment=isSecondaryAlignment,
                                           isNotPassingQualityControls=isNotPassingQualityControls,
                                           isUnmappedQuery=isUnmappedQuery),
                        which = GRanges(paste0(CHRSTRING, chr),
                                        IRanges(starts,ends)),
                        mapqFilter=mapqFilter)
    coverageTrack <- countBam(bamPath, param = sbp)
    coverageTrack$records <- coverageTrack$records
    return(coverageTrack)
}
