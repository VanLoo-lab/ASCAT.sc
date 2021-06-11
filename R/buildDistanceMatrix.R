buildDistanceMatrix <- function (meansSeg,
                                 weights,
                                 purs,
                                 ploidies,
                                 maxTumourPhi,
                                 distance="mse",
                                 gamma=1,
                                 sds=NULL,
                                 isPON=F,
                                 ismale=F,
                                 isX=F)
{
    errs <- matrix(NA, length(purs), length(ploidies))
    rownames(errs) <- purs
    colnames(errs) <- ploidies
    for (pp in 1:length(purs)) {
        for (pl in 1:length(ploidies)) {
            if (getTumourPhi(ploidies[pl], purs[pp]) > maxTumourPhi)
                errs[pp, pl] <- Inf
            else errs[pp, pl] <- geterrors(rho = purs[pp],
                                           phi = ploidies[pl],
                                           meansSeg=meansSeg,
                                           weights=weights,
                                           sds=if(distance=="mse") (pl/pp) * 0 + 1  else sds,
                                           gamma=gamma,
                                           distance=distance,
                                           isPON=isPON,
                                           ismale=ismale,
                                           isX=isX)
        }
    }
    errs
}

if(F)
    buildDistanceMatrix <-
        function(meansSeg,
                 weights,
                 purs,
                 ploidies,
                 maxTumourPhi)
    {
        step <- abs(ploidies[2]-ploidies[1])
        minp <- min(2,min(ploidies))
        maxp <- max(2,max(ploidies))
        averagepl <- seq(minp,maxp,step)
        errs <- matrix(NA, length(purs), length(averagepl))
        rownames(errs) <- purs
        colnames(errs) <- averagepl
        for (pp in 1:length(purs)) {
            for (pl in 1:length(averagepl)) {
                tumor_ploidy <- getTumourPhi(averagepl[pl],purs[pp])
                if(tumor_ploidy>maxTumourPhi | !any(abs(tumor_ploidy-ploidies)<step*2))
                    errs[pp,pl] <- Inf
                else
                    errs[pp, pl] <- geterrors(rho = purs[pp], phi = averagepl[pl],
                                              meansSeg, weights, (pl/pp) * 0 + 1)
            }
        }
        errs
    }
