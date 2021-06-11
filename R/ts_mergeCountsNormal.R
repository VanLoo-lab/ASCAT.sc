ts_mergeCountsNormal <- function(bamsN,
                                 window,
                                 lExclude,
                                 lSe = NULL,
                                 lGCT = NULL,
                                 allchr = 1:22,
                                 sdNormalise = 0,
                                 mc.cores=7,
                                 verbose=2)
{
    require(parallel)
    if(verbose>0) print("Get coverage in windows")
    lCT <- mclapply(bamsN,function(bamfile)
    {
        if(verbose>1) cat(".")
        lapply(allchr, function(chr) getCoverageTrack(bamPath = bamfile,
                                                      chr = paste0(chr),
                                                      lSe[[chr]]$starts,
                                                      lSe[[chr]]$ends))
    },mc.cores=mc.cores)
    if(verbose>0) print("Get coverage in targets")
    lCTex <- mclapply(bamsN,function(bamfile)
    {
        if(verbose>1)  cat(".")
        lapply(allchr, function(chr) getCoverageTrack(bamPath = bamfile,
                                                      chr = paste0(chr),
                                                      lExclude[[chr]]$starts,
                                                      lExclude[[chr]]$ends))
    },mc.cores=mc.cores)
    if(verbose>0) print("Remove target coverage from windows")
    lCTex <- lapply(1:length(lCT), function(x) ts_removeOnTargets(lCT[[x]],lCTex[[x]]))
    print("Add window coverage across BAMs")
    lCT. <- mclapply(1:length(allchr),function(chr)
    {
        if(verbose>1) cat(".")
        mm <- lCTex[[1]][[chr]]
        mm <- cbind(mm,lCTex[[1]][[chr]][,"records"])
        colnames(mm)[ncol(mm)] <- paste0("records_",bamsN[1])
        for(i in 2:length(lCT))
        {
            mm[,"records"] <- mm[,"records"]+lCTex[[i]][[chr]][,"records"]
            mm <- cbind(mm,lCTex[[i]][[chr]][,"records"])
            colnames(mm)[ncol(mm)] <- paste0("records_",bamsN[i])
        }
        mm
    },mc.cores=mc.cores)
    if(verbose>0) print("Return results")
    lCT.
}
