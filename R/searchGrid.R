searchGrid <- function (tracksSingle,
                        purs = seq(0.05, 1, 0.01),
                        ploidies = seq(1.7,
                                       6, 0.02),
                        forcepurity=NULL,
                        maxTumourPhi = 8,
                        errs = NULL,
                        radius=10,
                        localattempts=15,
                        maxattempts=150,
                        distance=c("mse", "statistical"),
                        ismale=F,
                        isPON=F,
                        gamma=1)
{
    errs. <- errs
    tracksSingle <- normaliseByPloidy(tracksSingle)
    meansSeg <- unlist(lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- unlist(lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
                return(NA)
            mu <- mean(tracksSingle$lCTS[[i]]$smoothed[isIn],
                       na.rm = T)
            mu
        }))
    }))
    sds <- NULL
    isX <- FALSE
    if(ismale & isPON)
        isX <- unlist(lapply(1:length(tracksSingle$lSegs), function(i) {
            out <- tracksSingle$lSegs[[i]]$output
            isx <- unlist(lapply(1:nrow(out), function(x) {
                if(out$chrom%in%c("X","23")) return(TRUE)
                return(FALSE)
            }))
        }))
    if(grepl("stat",distance[1]))
        sds <- unlist(lapply(1:length(tracksSingle$lSegs), function(i) {
            out <- tracksSingle$lSegs[[i]]$output
            means <- unlist(lapply(1:nrow(out), function(x) {
                isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                    tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
                if (sum(isIn) < 2)
                    return(NA)
                sds <- sd(tracksSingle$lCTS[[i]]$smoothed[isIn],
                          na.rm = T)
                sds
            }))
        }))
    weights <- unlist(lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output$num.mark
    }))
    if (is.null(errs))
        errs <- buildDistanceMatrix(meansSeg,
                                    weights,
                                    purs,
                                    ploidies,
                                    maxTumourPhi,
                                    gamma=gamma,
                                    distance=distance[1],
                                    ismale=ismale,
                                    isPON=isPON,
                                    isX=isX,
                                    sds=sds)
    SOL_SEARCHING <- T
    counts <- 0
    totalcounts <- 0
    ambiguous <- F
    while (SOL_SEARCHING)
    {
        mins <- arrayInd(which.min(errs), dim(errs))
        purity <- purs[mins[1]]
        ploidy <- ploidies[mins[2]]
        meansSeg <- fitProfile(tracksSingle,
                               purity,
                               ploidy,
                               gamma=gamma,
                               ismale=ismale,
                               isPON=isPON)
        SOL_SEARCHING <- !isGoodSolution(meansSeg)
        if (SOL_SEARCHING)
        {
            errs[mins[1],mins[2]] <-  Inf
            counts <- counts+1
            totalcounts <- totalcounts+1
            if(counts>localattempts)
            {
                minX <- max(mins[1]-radius,1)
                maxX <- min(mins[1]+radius,nrow(errs))
                minY <- max(mins[2]-radius,1)
                maxY <- min(mins[2]+radius,ncol(errs))
                errs[minX:maxX,minY:maxY] <- Inf
                counts <- 0
            }
            if(totalcounts>maxattempts)
            {
                SOL_SEARCHING <- FALSE
                ambiguous <- T
            }
        }
    }
    cat(if(ambiguous) "Too many local attempts - solution not found and sample/cell marked as ambiguous\n" else "Solution found.\n")
    if(is.null(forcepurity))
        return(list(errs = errs,
                    purity = purity,
                    ploidy = ploidy,
                    ambiguous = ambiguous))
    else
    {
        if(all(abs(purity-forcepurity)<=0.02))
            return(list(errs = errs,
                        purity = purity,
                        ploidy = ploidy,
                        ambiguous = ambiguous))
        else
        {
            return(append(searchGrid(tracksSingle=tracksSingle,
                                     purs = forcepurity,
                                     ploidies = ploidies,
                                     forcepurity = NULL,
                                     maxTumourPhi = maxTumourPhi,
                                     errs = errs.,
                                     radius = radius,
                                     localattempts = localattempts,
                                     maxattempts = maxattempts,
                                     gamma=gamma,
                                     ismale=ismale,
                                     isPON=isPON,
                                     distance=distance),
                          list(bestfit=list(errs = errs,
                                            purity = purity,
                                            ploidy = ploidy,
                                            ambiguous = ambiguous))))
        }
    }
}
