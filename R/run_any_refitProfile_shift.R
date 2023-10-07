run_any_refitProfile_shift <- function(res,
                                       sample_indice,
                                       shift=1,
                                       CHRS=c(1:22,"X","Y"),
                                       outdir,
                                       print=T,
                                       gridpur=seq(-.05,.05,.01),
                                       gridpl=seq(-.1,.2,.01))
{
    if(!"allProfiles.refitted.manual"%in%names(res))
    {
        res$allProfiles.refitted.manual <- list()
    }
    if(is.null(res$isPON)) res$isPON <- F
    gamma <- if("gamma"%in%names(res)) res$gamma else if("GAMMA"%in%names(res)) res$GAMMA else 1
    solution <- res$allSolutions.refitted.manual[[sample_indice]] <- refitProfile_shift(track=res$allTracks.processed[[sample_indice]],
                                                                                        solution=res$allSolutions.refitted.auto[[sample_indice]],
                                                                                        gamma=gamma,
                                                                                        ismale=res$sex[sample_indice]=="male",
                                                                                        isPON=res$isPON,
                                                                                        CHRS=res$chr,
                                                                                        shift=shift,
                                                                                        gridpur=gridpur,
                                                                                        gridpl=gridpl)
    res$allProfiles.refitted.manual[[sample_indice]] <- getProfile(fitProfile(res$allTracks.processed[[sample_indice]],
                                                                              solution$purity,
                                                                              solution$ploidy,
                                                                              gamma=gamma,
                                                                              ismale=res$sex[sample_indice]=="male",
                                                                              isPON=res$isPON),
                                                                   CHRS=res$chr)
    writeProfile(res$allProfiles.refitted.manual[[sample_indice]],
                 samplename=paste0(names(res$allTracks.processed)[sample_indice],"-manual_refit"),
                 outdir=outdir)
    res
}
