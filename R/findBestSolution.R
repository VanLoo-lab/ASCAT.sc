## ##################################################
## Deprecated
## ##################################################
findBestSolution <-function (tracksSingle,
                             colSeg = rgb(0.5, 0.2, 0.5, 0.7),
                             lwdSeg = 2,
                             REFs = c(1:23),
                             purs=seq(0.05, 1, 0.005),
                             ploidies=seq(1.7, 6, 0.01),
                             maxTPhi=5,
                             gamma=1,
                             ...)
{
    meanS <- mean(unlist(lapply(REFs, function(x) {
        oo <- tracksSingle$lSegs[[x]]$output
        log2(10) * inverse.rle(list(values = oo$seg.mean, lengths = oo$num.mark))
    })), na.rm = T)
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - meanS
        return(x)
    })
    tracksSingle$lSeg <- lapply(tracksSingle$lSeg, function(x) {
        x$output$seg.mean <- log2(10) * x$output$seg.mean - meanS
        return(x)
    })
    searchGrid(tracksSingle,
               purs,
               ploidies,
               maxTumourPhi,
               gamma=gamma)
}
