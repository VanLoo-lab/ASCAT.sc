getCoverageTrack.Fix <-
function(bamPath,
                                 chr,
                                 lengthChr,
                                 step,
                                 CHRSTRING="")
{
    require(Rsamtools)
    require(GenomicRanges)
    divideChr <- seq(0, lengthChr, step)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE),
        which = GRanges(paste0(CHRSTRING,chr), IRanges(starts, ends)))
    coverageTrack <- countBam(bamPath, param = sbp)
    return(coverageTrack)
  }
