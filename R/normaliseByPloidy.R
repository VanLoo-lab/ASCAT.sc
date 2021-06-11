normaliseByPloidy <- function(tracksSingle)
{
    meanS <- mean(unlist(lapply(1:length(tracksSingle$lSegs),
                                function(x) {
                                    oo <- tracksSingle$lSegs[[x]]$output
                                    inverse.rle(list(values = oo$seg.mean, lengths = oo$num.mark))
                                })), na.rm = T)
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- x$smoothed - meanS
        return(x)
    })
    tracksSingle$lSeg <- lapply(tracksSingle$lSeg, function(x) {
        x$output$seg.mean <- x$output$seg.mean - meanS
        return(x)
    })
    tracksSingle
}
