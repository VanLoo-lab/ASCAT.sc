getTrackForAll.bins <- function (logr, CHRS, STARTS, ENDS, allchr=ALLCHR)
{
    lT <- lapply(allchr, function(chr)
    {
        cond <- CHRS==paste0("chr",chr)
        data.frame(start=STARTS[cond],
                   end=ENDS[cond],
                   records=logr[cond],
                   fitted=logr[cond],
                   smoothed=logr[cond])
    })
    lSe <- lapply(allchr, function(chr)
    {
        cond <- CHRS==paste0("chr",chr)
        list(starts=STARTS[cond],ends=ENDS[cond])
    })
    lSegs <- lapply(1:length(lT), function(x) {
        require(DNAcopy)
        segments <- segmentTrack(lT[[x]]$smoothed,
                                 chr = paste0(x),
                                 sd = 0,
                                 lSe[[x]]$starts,
                                 lSe[[x]]$ends,
                                 ALPHA=0.001,
                                 smooth=F)
    })
    names(lSegs) <- paste0(1:length(lT))
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
}
