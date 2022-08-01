getLSegs.multipcf <- function(allTracks,
                              lCTS,
                              lSe,
                              lGCT,
                              lNormals=NULL,
                              allchr,
                              segmentation_alpha=0.01,
                              MC.CORES=1)
{
    ## #############################################################
    require(parallel)
    ## #############################################################
    smoothAll <- function(lCTS, lSe, lGCT, lNormals, allchr, MC.CORES)
    {
        lCTSs <- parallel::mclapply(lCTS, function(lCT)
        {
            cat(".")
            lCTS <- smoothCoverageTrack(lCT=lCT,lSe=lSe,lGCT=lGCT, lNormals=lNormals)
            names(lCTS) <- allchr
            lCTS
        },mc.cores=MC.CORES)
        cat("\n")
        lCTSs
    }
    ## #############################################################
    print("Smoothing all tracks")
    lCTSs <- smoothAll(lCTS,lSe, lGCT, lNormals, allchr, MC.CORES=MC.CORES)
    ## #############################################################
    runMultiPCF <- function(allT, penalties, nchr, mc.cores=1)
    {
        require(copynumber)
        chr_pcfed <- mclapply(1:nchr,function(i)
        {
            cat(".")
            tmpdata <- data.frame("chr"=rep(i,nrow(allT$lCTS[[1]][[i]])),
                                  pos=round(allT$lCTS[[1]][[i]][,"start"]+allT$lCTS[[1]][[i]][,"width"]/2),
                                  do.call("cbind",lapply(allT$lCTS,function(x)
                                  {
                                      x[[i]]$smoothed
                                  })))
            out <- lapply(penalties,function(penalty)
            {
                invisible(capture.output(suppressMessages({ok <- copynumber::multipcf(data=tmpdata, gamma=penalty)})))
                ok
            })
            out
        },mc.cores=mc.cores)
        names(chr_pcfed) <- paste0("chr",1:nchr)
        chr_pcfed
    }
    ## #############################################################
    nchr <- length(lCTS[[1]])
    allT <- list(lCTS=lCTSs)
    penalties <- 1/segmentation_alpha ## rough translation from CBS pval to pcf penalty
    ## #############################################################
    print("Segmenting all chromosomes across tracks")
    allsegs <- runMultiPCF(allT,
                           penalties=penalties,
                           nchr=nchr,
                           mc.cores=MC.CORES)
    ## #############################################################
    allTracks.processed <- mclapply(1:length(allTracks), function(i)
    {
        allTracks[[i]]$lSegs <- lapply(allsegs,function(x)
        {
            x <- x[[1]]
            list(data=NULL,
                 output=data.frame(ID="Sample.ID",
                                   chrom=x[,"chrom"],
                                   loc.start=as.numeric(as.character(x[,"start.pos"])),
                                   loc.end=as.numeric(as.character(x[,"end.pos"])),
                                   num.mark=as.numeric(as.character(x[,"n.probes"])),
                                   seg.mean=as.numeric(as.character(x[,5+i]))))
        })
        names(allTracks[[i]]$lSegs) <- allchr
        allTracks[[i]]$lCTS <- lCTSs[[i]]
        names(allTracks[[i]]$lCTS) <- allchr
        allTracks[[i]]
    },mc.cores=MC.CORES)
    names(allTracks.processed) <- names(allTracks)
    allTracks.processed
}
