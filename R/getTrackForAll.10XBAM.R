getTrackForAll.10XBAM <- function (bamfile,
                                   window,
                                   REFGENOME,
                                   lSe = NULL,
                                   lGCT = NULL,
                                   allchr = 1:22,
                                   CHRSTRING="chr",
                                   sdNormalise = 0,
                                   segmentation_alpha = 0.01,
                                   isDuplicate = F,
                                   isSecondaryAlignment = F,
                                   isNotPassingQualityControls = NA,
                                   isUnmappedQuery = NA,
                                   mapqFilter = 0,
                                   barcodes=NULL)
{
    if (is.null(lSe)) {
        print("get Start-End of segments")
        lSe <- lapply(allchr, function(chr) getStartsEnds(window,
                                                          paste0(chr)))
    }
    if(is.null(barcodes))
    {
        print("get Barcodes")
        barcodes <- getCoverageTrack.10XBAM(bamPath = bamfile,
                                            chr = allchr[1],
                                            lSe[[allchr[1]]]$starts,
                                            lSe[[allchr[1]]]$ends,
                                            CHRSTRING=CHRSTRING,
                                            isDuplicate = isDuplicate,
                                            isSecondaryAlignment = isSecondaryAlignment,
                                            isNotPassingQualityControls = isNotPassingQualityControls,
                                            isUnmappedQuery = isUnmappedQuery, mapqFilter = mapqFilter,
                                            barcodes=NULL)$barcodes
    }
    print("get Coverage Track")
    lCT <- lapply(allchr, function(chr)
        getCoverageTrack.10XBAM(bamPath = bamfile,
                                chr = paste0(chr),
                                lSe[[chr]]$starts,
                                lSe[[chr]]$ends,
                                CHRSTRING=CHRSTRING,
                                isDuplicate = isDuplicate,
                                barcodes=barcodes,
                                isSecondaryAlignment = isSecondaryAlignment,
                                isNotPassingQualityControls = isNotPassingQualityControls,
                                isUnmappedQuery = isUnmappedQuery, mapqFilter = mapqFilter))
    if (is.null(lGCT)) {
        print("get GC content")
        lGCT <- lapply(allchr, function(chr) gcTrack(chr, lSe[[chr]]$starts,
                                                     lSe[[chr]]$ends, dna=REFGENOME))
    }
    lCTs <- lapply(1:length(lCT[[1]][[1]]), function(cell)
        lapply(1:length(allchr),
               function(chr) lCT[[chr]][[1]][[cell]]))
    print("correct for GC content")
    lCTSs <- lapply(lCTs, function(lCT) smoothCoverageTrack(lCT, lSe, lGCT))
    gc()
    print("segment Tracks")
    lSegs <- lapply(lCTSs,function(lCTS)
    {
        lSegs <- lapply(1:length(lCTS), function(x)
        {
            require(DNAcopy)
            invisible(capture.output(segments <- segmentTrack(lCTS[[x]]$smoothed, chr = paste0(x),
                                                              sd = sdNormalise, lSe[[x]]$starts, lSe[[x]]$ends,
                                                              ALPHA = segmentation_alpha)))
            segments
        })
        names(lSegs) <- paste0(1:length(lCT))
        lSegs
    })
    alltracks  <- list(lCTS = lCTSs, lSegs = lSegs)
    return(alltracks)
}
