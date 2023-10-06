run_any_refitProfile <- function(res,
                                 sample_indice,
                                 chr1,
                                 ind1,
                                 n1,
                                 chr2,
                                 ind2,
                                 n2,
                                 CHRS=c(1:22,"X","Y"),
                                 gridpur=seq(-.05,.05,.01),
                                 gridpl=seq(-.1,.2,.01))
{
    if(!"allProfiles.refitted.manual"%in%names(res))
    {
        res$allProfiles.refitted.manual <- list()
    }
    if(is.null(res$isPON)) res$isPON <- F
    gamma <- if("gamma"%in%names(res)) res$gamma else if("GAMMA"%in%names(res)) res$GAMMA else 1
    solution <- res$allSolutions.refitted.manual[[sample_indice]] <- refitProfile(track=res$allTracks.processed[[sample_indice]],
                                                                 solution=res$allSolutions.refitted.auto[[sample_indice]],
                                                                 chr1=chr1,
                                                                 ind1=ind1,
                                                                 total1=n1,
                                                                 chr2=chr2,
                                                                 ind2=ind2,
                                                                 total2=n2,
                                                                 gamma=gamma,
                                                                 ismale=res$sex[sample_indice]=="male",
                                                                 isPON=res$isPON,
                                                                 CHRS=res$chr,
                                                                 gridpur=gridpur,
                                                                 gridpl=gridpl)
    res$allProfiles.refitted.manual[[sample_indice]] <- getProfile(fitProfile(res$allTracks.processed[[sample_indice]],
                                                                               solution$purity,
                                                                               solution$ploidy,
                                                                               gamma=gamma,
                                                                               ismale=res$sex[sample_indice]=="male",
                                                                               isPON=res$isPON),
                          CHRS=res$chr)
    res
}
