predictRefit_all <- function(res, ismedian=FALSE, gamma=1)
{
    preds <- lapply(res$allProfiles, function(x) try(predictRefit(x)))
    res$allProfiles.refitted.auto <- list()
    res$allSolutions.refitted.auto <- list()
    for(i in 1:length(preds))
    {
        if(is.numeric(preds[[i]]) & length(preds[[i]])==1)
        {
            if(preds[i]!=1)
                {
                    res$allSolutions.refitted.auto[[i]] <- try(refitProfile_shift(track=res$allTracks.processed[[i]],
                                                                                  solution=res$allSolutions[[i]],
                                                                                  CHRS=res$chr,
                                                                                  gamma=gamma,
                                                                                  ismedian=ismedian,
                                                                                  shift=if(preds[i]==0) 1 else -1,
                                                                                  isPON=res$isPON),silent=F)
                    res$allProfiles.refitted.auto[[i]] <- try(getProfile(fitProfile(res$allTracks.processed[[i]],
                                                                                    purity=res$allSolutions.refitted.auto[[i]]$purity,
                                                                                    ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                                                                    ismedian=ismedian,
                                                                                    gamma=gamma,
                                                                                    ismale=if(res$sex[i]=="male") T else F),
                                                                         CHRS=res$chr),silent=F)
                }
            if(preds[i]==1)
            {
                res$allSolutions.refitted.auto[[i]] <- res$allSolutions[[i]]
                res$allProfiles.refitted.auto[[i]] <- try(getProfile(fitProfile(res$allTracks.processed[[i]],
                                                                                purity=res$allSolutions.refitted.auto[[i]]$purity,
                                                                                ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                                                                gamma=gamma,
                                                                                ismedian=ismedian,
                                                                                ismale=if(res$sex[i]=="male") T else F),
                                                                     CHRS=res$chr),silent=F)
            }
        }
    }
    return(res)
}
