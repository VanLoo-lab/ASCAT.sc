getTrackForAll.bins <- function (logr,
                                 CHRS,
                                 STARTS,
                                 ENDS,
                                 allchr=ALLCHR.,
                                 segmentation_alpha=0.001,
                                 transform=FALSE,
                                 ismedian=FALSE,
                                 min.width=5,
                                 SBDRY=NULL)
{
    suppressPackageStartupMessages(require(DNAcopy))
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
        segments <- segmentTrack_pcf(lT[[x]]$smoothed,
                                     chr = paste0(x),
                                     sd = 0,
                                     lSe[[x]]$starts,
                                     lSe[[x]]$ends,
                                     transform=transform,
                                     ismedian=ismedian,
                                     smooth=FALSE,
                                     ALPHA=segmentation_alpha,
                                     min.width=min.width,
                                     SBDRY=SBDRY)
    })
    names(lSegs) <- names(lT) <- allchr
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
}
