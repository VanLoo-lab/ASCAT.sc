fitProfile <- function(tracksSingle,
                       purity,
                       ploidy,
                       ismale=F,
                       isPON=F,
                       gamma=1,
                       gridsearch=F)
{
    tracksSingle <- normaliseByPloidy(tracksSingle)
    NN <- length(tracksSingle$lSegs)
    if(gridsearch) NN <- min(NN,23) ## avoid checking chromosome Y in the grid search (leads to HD flag in males)
    meansSeg <- lapply(1:NN, function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start >= out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 1)
                return(list(roundmu = NA, mu = NA, sd = NA, start = out$loc.start[x],
                            end = out$loc.end[x],
                            num.mark=out$num.mark[x],
                            num.mark.in=sum(isIn)))
            mu <- mean(tracksSingle$lCTS[[i]]$smoothed[isIn],
                         na.rm = T)
            if (sum(isIn) < 2)
                sd <- NA
            else
                sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                          na.rm = T)
            list(roundmu = transform_bulk2tumour(mu,
                                                 purity,
                                                 ploidy,
                                                 gamma=gamma,
                                                 ismale=ismale,
                                                 isPON=isPON,
                                                 isX=tracksSingle$lSegs[[i]]$output$chrom[1]%in%c("23","X", "chrX","chr23")),
                 mu = mu,
                 sd = sd,
                 start = out$loc.start[x],
                 end = out$loc.end[x],
                 num.mark=out$num.mark[x],
                 num.mark.in=sum(isIn))
        })
    })
}
