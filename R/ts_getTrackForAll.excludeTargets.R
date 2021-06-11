ts_getTrackForAll.excludeTargets <- function (bamfile,
                                              window,
                                              lExclude,
                                              lCTS.normal,
                                              lSe,
                                              lGCT,
                                              lRepli = NULL,
                                              allchr = 1:22,
                                              sdNormalise = 0)
{
    print("get Coverage Track")
    lCT <- lapply(allchr, function(chr) getCoverageTrack(bamPath = bamfile,
                                                         chr = paste0(chr), lSe[[chr]]$starts, lSe[[chr]]$ends))
    print("get Coverage Track in Excluded regions")
    lCTex <- lapply(allchr, function(chr) getCoverageTrack(bamPath = bamfile,
                                                           chr = paste0(chr), lExclude[[chr]]$starts, lExclude[[chr]]$ends))
    lCT <- ts_removeOnTargets(lCT,lCTex)
    print("correct for GC, replication timing and other biases in normal")
    lCTS <- ts_smoothCoverageTrackAll(lCT,lSe=lSe, lGCT=lGCT, lRepli=lRepli,lNormal=lCTS.normal)
    gc()
    print("segment Tracks")
    lSegs <- lapply(1:length(lCTS), function(x) {
        require(DNAcopy)
        segments <- segmentTrack(lCTS[[x]]$smoothed, chr = paste0(x),
                                 sd = sdNormalise, lSe[[x]]$starts, lSe[[x]]$ends)
    })
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS = lCTS, lSegs = lSegs)
    return(tracks)
}
