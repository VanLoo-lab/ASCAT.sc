get_filters_targeted_sequencing <- function(normal_bams,
                                            fasta,
                                            ismale=c(T,T,T,F,F,F),
                                            allchr=c(1:22,"X","Y", "MT"),
                                            autosomes=1:22,
                                            chrX=23,
                                            chrY=24,
                                            chrMT=25,
                                            pcremove=15,
                                            pcremoveY=60,
                                            pcremoveMT=15,
                                            MC.CORES=1)
{
    require(parallel)
    require(Rsamtools)
    require(Biostrings)
    START_WINDOW <- 10000
    print("## load genome track")
    dna <- getRefGenome(fasta=fasta, CHRS=allchr)
    print("## get bins ")
    lSe <- mclapply(ALLCHR,function(chr)
    {
        getStartsEnds(window=START_WINDOW,
                      chr=paste0("",chr),
                      exclude=bed,
                      dna=dna,
                      lengthChr=sapply(dna,length)[chr])
    }, mc.cores=MC.CORES)
    names(lSe) <- ALLCHR
    print("## get GC in bins ")
    lGCT <- lapply(allchr,function(chr)
    {
        gcT <- gcTrack(chr,
                       lSe[[chr]]$starts,
                       lSe[[chr]]$ends,
                       dna=dna)
    })
    names(lSe) <- names(lGCT) <- allchr
    print("## get GC in bins ")
    allTracks <- lapply(normal_bams,function(x)
    {
        try(getTrackForAll(x,
                           START_WINDOW,
                           lSe=lSe,
                           lGCT=lGCT,
                           allchr=ALLCHR,
                           mapqFilter=30,
                           sdNormalise=0),silent=T)
    },mc.cores=MC.CORES)
    getfilters <- function(allTracks,
                           sex,
                           autosomes=1:22,
                           chrX=23,
                           chrY=24,
                           chrMT=25,
                           pcremove=15,
                           pcremoveY=60,
                           pcremoveMT=15)
    {
        getremove <- function(mmat, pc)
        {
            ranks <- apply(mmat, 2, function(x) rank(x))
            NN <- nrow(ranks)
            sums <- rowSums(apply(ranks,2,function(x) log(x/(NN+1-x))))
            ord <- order(abs(sums),decreasing=T)
            remove <- ord[1:round(length(ord)*pc/100)]
        }
        mmat <- sapply(allTracks,function(x)
        {
            do.call("c",lapply(autosomes,function(y) x$lCTS[[y]]$smoothed))
        })
        ## ##################################################
        removeAutosome <- getremove(mmat, pc=pcremove)
        ## ##################################################
        mmat <- sapply(1:length(allTracks),function(x)
        {
            do.call("c",lapply(chrX,function(y) allTracks[[x]]$lCTS[[y]]$smoothed+ifelse(ismale,1,0)))
        })
        ## ##################################################
        removeX <- getremove(mmat, pc=pcremove)
        ## ##################################################
        mmat <- sapply(1:length(allTracks),function(x)
        {
            do.call("c",lapply(chrY,function(y) allTracks[[x]]$lCTS[[y]]$smoothed+ifelse(ismale,1,2)))
        })
        ## ##################################################
        removeY <- getremove(mmat, pc=pcremoveY)
        ## ##################################################
        mmat <- sapply(1:length(allTracks),function(x)
        {
            do.call("c",lapply(chrMT,function(y) allTracks[[x]]$lCTS[[y]]$smoothed))
        })
        mmat <- t(t(mmat)-apply(mmat,2,median, na.rm=T))
        ## ##################################################
        removeMT <- getremove(mmat, pc=pcremoveMT)
        ## ##################################################
        indices <- do.call("rbind",
                           lapply(1:length(lSe), function(i)
                               do.call("rbind",lapply(1:length(lSe[[i]]$start),function(j) c(i,j)))))
        ## ##################################################
        indices_remove <- rbind(indices[remove,],
                                indices[indices[,1]==chrX,][removeX,],
                                indices[indices[,1]==chrY,][removeY,],
                                indices[indices[,1]==chrMT,][removeMT,])
        ## ##################################################
        nlSe <- lapply(1:length(lSe),function(x)
        {
            keep <- !1:length(lSe[[x]]$starts)%in%indices_remove[indices_remove[,1]==x,2]
            list(starts=lSe[[x]]$starts[keep],
                 ends=lSe[[x]]$ends[keep])
        })
        ## ##################################################
        filters <- list(lSe=lSe,
                        lSe.filtered=nlSe,
                        indices_remove=indices_remove)
    }
    getfilters(allTracks=allTracks,
               sex=sex,
               autosomes=autosomes,
               chrX=chrX,
               chrY=chrxY,
               chrMT=chrMT,
               pcremove=pcremove,
               pcremoveY=pcremoveY,
               pcremoveMT=pcremoveMT)
}
