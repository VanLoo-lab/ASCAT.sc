predictRefit_all <- function(res, ismedian=FALSE, gamma=1)
{
    require(xgboost)
    model_path <- system.file("extdata",
                              "model_xgboost_tcga_snp6_200426.json",
                              package = "ASCAT.sc")
    if(model_path=="")
    {
        stop("Model file not found in package. Check installation.")
    }
    xgb_model <- xgboost::xgb.load(model_path)
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
