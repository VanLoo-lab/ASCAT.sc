getCoverageTrack.10XBAM <- function (bamPath,
                                     chr,
                                     starts,
                                     ends,
                                     CHRSTRING = "",
                                     isDuplicate=F,
                                     isSecondaryAlignment=F,
                                     isNotPassingQualityControls=NA,
                                     isUnmappedQuery=NA,
                                     mapqFilter=0,
                                     barcodes=NULL,
                                     tag=c("CB","CR"))
{
    require(Rsamtools)
    require(GenomicRanges)
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate=isDuplicate,
                                           isSecondaryAlignment=isSecondaryAlignment,
                                           isNotPassingQualityControls=isNotPassingQualityControls,
                                           isUnmappedQuery=isUnmappedQuery),
                        which = GRanges(paste0(CHRSTRING, chr),
                                        IRanges(starts,ends)),
                        mapqFilter=mapqFilter,
                        tag=tag[1])
    coverageTrack <- scanBam(bamPath, param = sbp)
    .guessBarcodes <- function(cT, pcchromosome=.8)
    {
        counts <- table(unlist(lapply(cT,function(x) unique(unlist(x)))))
        names(counts)[counts>length(cT)*pcchromosome]
    }
    if(is.null(barcodes)) barcodes <- .guessBarcodes(coverageTrack)
    .countBarcodes <- function(barcodesBAM, barcodesUNIQUE)
    {
        counts <- as.vector(table(unlist(barcodesBAM))[barcodesUNIQUE])
        names(counts) <- barcodesUNIQUE
        counts[is.na(counts)] <- 0
        counts
    }
    .formatCounts <- function(counts,chr,start,end)
    {
        data.frame(space=NA,
                   start=start,
                   end=end,
                   width=NA,
                   file=NA,
                   records=counts,
                   nucleotides=NA)
    }
    .formatOutput <- function(cT, chr, starts, ends)
    {
        out <- lapply(cT,.countBarcodes,barcodes)
        out <- lapply(barcodes,
                      function(bc)
                      {
                          df <- do.call("rbind",
                                        lapply(1:length(out),
                                               function(i)
                                                   .formatCounts(out[[i]][bc],chr,starts[i],ends[i])))
                          df$width <- df$end-df$start
                          df$space <- chr
                          df
                      })
        names(out) <- barcodes
        out
    }
    coverageTrack <- .formatOutput(coverageTrack,chr,starts,ends)
    return(list(coverageTrack=coverageTrack,
                barcodes=barcodes))
}
