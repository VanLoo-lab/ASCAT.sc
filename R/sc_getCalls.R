sc_getCalls <- function (tracksSingle, purity, ploidy)
{
    tracksSingle <- normaliseByPloidy(tracksSingle)
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
                return(list(roundmu = NA, mu = NA, sd = NA, start = out$loc.start[x],
                            end = out$loc.end[x]))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
                         na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                      na.rm = T)
            list(roundmu = transform_bulk2tumour(mu, purity, ploidy), mu = mu,
                 sd = sd, start = out$loc.start[x], end = out$loc.end[x])
        })
        t <- cbind(out,
                   sapply(means, function(y) y$roundmu),
                   sapply(means, function(y) y$mu),
                   sapply(means, function(y) y$sd))
        colnames(t) <- c(colnames(out), "CNtot", "CNseg", "CNsd")
        t
    })
    return(meansSeg)
}
