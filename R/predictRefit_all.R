predictRefit_all <- function(res, ismedian=FALSE, gamma=1)
{
    if(is.null(res$is_methyl))
    {
        data("model_xgboost",package="ASCAT.sc")
        xgb_model <- final_model
    }
    else
    {
        data("model_xgb_trained_on_tcga_toby",package="ASCAT.ma")
        xgb_model <- model_xgb
    }
    predict_path <-function(df)
    {
        fit_df <- get_fit_features(df)
        good_fit_probability <- predict(xgb_model, newdata = as.matrix(fit_df))
        return(good_fit_probability)
    }
    get_optimal_solution <-function(sol, res, x, N=10)
    {
        errs <- sol$errs
        bao<-findLocalMinima(errs, N=N)$bao
        prof<-list()
        for(i in 1:N)
        {
                purity<-as.numeric(rownames(errs)[bao[i,1]])
                ploidy<-as.numeric(colnames(errs)[bao[i,2]])
                prof[[i]]<-getProfile(fitProfile(res$allTracks.processed[[x]],
                                                 purity=purity,
                                                 ploidy=ploidy,
                                                 ismedian=ismedian,
                                                 gamma=gamma,
                                                 ismale=if(res$sex[x]=="male") T else F),
                                      CHRS=res$chr[!res$chr%in%c("X","Y","chrX","chrY")])
        }
        prds<-sapply(prof,function(x)
        {
            predict_path(x)
        })
        best_i <- which.min(prds)[1]
        nsol <- sol
        nsol$purity <- as.numeric(rownames(errs)[bao[best_i,1]])
        nsol$ploidy <- as.numeric(colnames(errs)[bao[best_i,2]])
        nsol
    }

    for(i in 1:length(res$allSolutions))
    {
        res$allSolutions.refitted.auto[[i]] <- try(get_optimal_solution(res$allSolutions[[i]],
                                                                            res,
                                                                            x=i),silent=F)
        res$allProfiles.refitted.auto[[i]] <- try(getProfile(fitProfile(res$allTracks.processed[[i]],
                                                                        purity=res$allSolutions.refitted.auto[[i]]$purity,
                                                                        ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                                                        ismedian=ismedian,
                                                                        gamma=gamma,
                                                                        ismale=if(res$sex[i]=="male") T else F),
                                                             CHRS=res$chr),silent=F)
    }
    return(res)
}

if(FALSE){
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
}
