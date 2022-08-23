getTrackForAll.10XBAM <- function (bamfile,
                                   window=NULL,
                                   refgenome=NULL,
                                   lSe = NULL,
                                   lGCT = NULL,
                                   allchr = 1:22,
                                   chrstring="chr",
                                   sdNormalise = 0,
                                   doSeg=TRUE,
                                   segmentation_alpha = 0.01,
                                   isDuplicate = F,
                                   isSecondaryAlignment = F,
                                   isNotPassingQualityControls = NA,
                                   isUnmappedQuery = NA,
                                   mapqFilter = 30,
                                   barcodes=NULL)
{
    if (is.null(lSe)) {
        if(is.null(window)) stop("either provide window or starts/ends")
        print("get Start-End of segments")
        lSe <- lapply(allchr, function(chr) getStartsEnds(window,
                                                          paste0(chr)))
    }
    if(is.null(barcodes) | is.na(barcodes))
    {
        print("get Barcodes")
        barcodes <- getCoverageTrack.10XBAM(bamPath = bamfile,
                                            chr = allchr[1],
                                            lSe[[allchr[1]]]$starts,
                                            lSe[[allchr[1]]]$ends,
                                            chrstring=chrstring,
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
                                chrstring=chrstring,
                                isDuplicate = isDuplicate,
                                barcodes=barcodes,
                                isSecondaryAlignment = isSecondaryAlignment,
                                isNotPassingQualityControls = isNotPassingQualityControls,
                                isUnmappedQuery = isUnmappedQuery, mapqFilter = mapqFilter))
    lCTs <- lapply(1:length(lCT[[1]][[1]]), function(cell)
    {
        lCTS <- lapply(1:length(allchr),
                    function(chr) lCT[[chr]][[1]][[cell]])
        names(lCTS) <- allchr
        lCTS
    })
    names(lCTs) <- barcodes
    if(!doSeg)
        return(lCTs)
    if (is.null(lGCT))
    {
        if(is.null(refgenome)) stop("provide either reference genome or GC content track")
        print("get GC content")
        lGCT <- lapply(allchr, function(chr) gcTrack(chr, lSe[[chr]]$starts,
                                                     lSe[[chr]]$ends, dna=refgenome))
    }
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
