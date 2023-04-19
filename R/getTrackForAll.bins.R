getTrackForAll.bins <- function (logr,
                                 CHRS,
                                 STARTS,
                                 ENDS,
                                 allchr=ALLCHR.,
                                 segmentation_alpha=0.001,
                                 min.width=5,
                                 SBDRY=NULL)
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
                                 transform=TRUE,
                                 ALPHA=segmentation_alpha,
                                 min.width=min.width,
                                 SBDRY=SBDRY)
    })
    names(lSegs) <- names(lT) <- allchr
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
}
