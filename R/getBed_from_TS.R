getBed_from_TS <- function(normal_bams,
                           build="hg19",
                           lSe=NULL,
                           output_filename="mybed_manual.bed",
                           allchr = c(1:22,"X","Y"),
                           MC.CORES=1)
{
    res <- list()
    print("## load bin starts and ends")
    if(build=="hg19")
    {
        data("lSe_filtered_30000.hg19",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg19",package="ASCAT.sc")
        allchr. <- gsub("chr","",allchr)
        res$lSe <- lapply(allchr., function(chr) lSe.hg19.filtered[[chr]])
        names(res$lSe) <- allchr
        names(lGCT.hg19.filtered) <- names(res$lSe)
        res$lGCT <- lapply(allchr, function(chr) lGCT.hg19.filtered[[chr]])
    }
    if(build=="hg38")
    {
        data("lSe_filtered_30000.hg38",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg38",package="ASCAT.sc")
        names(lGCT.hg38.filtered) <- names(lSe.hg38.filtered)
        allchr. <- paste0("chr",gsub(chrstring_bam,"",allchr))
        res$lSe <- lapply(allchr., function(chr) lSe.hg38.filtered[[chr]])
        names(res$lSe) <- allchr
        res$lGCT <- lapply(allchr., function(chr) lGCT.hg38.filtered[[chr]])
        names(res$lGCT) <- allchr
        if(chrstring_bam=="")
            names(res$lGCT) <- names(res$lSe) <- gsub("chr","",names(res$lSe))
    }
    if (build == "mm39")
    {
        data("lSe_unfiltered_5000.mm39", package = "ASCAT.sc")
        data("lGCT_unfiltered_5000.mm39", package = "ASCAT.sc")
        names(lGCT)[1:length(allchr)] <- names(lSe)[1:length(allchr)] <- allchr
        res$lSe <- lapply(allchr, function(x) lSe[[x]])
        res$lGCT <- lapply(allchr, function(x) lGCT[[x]])
        names(res$lGCT) <- names(res$lSe) <- allchr
    }
    if (!build %in% c("hg19", "hg38","mm39"))
    {
        START_WINDOW <- median(unlist(lapply(lSe,function(x) x$ends-x$starts+1)))
        res$lSe <- lSe
        res$lGCT <- lGCT
        names(res$lGCT) <- names(res$lSe) <- allchr
    }
    lSe <- res$lSe
    split_lSe <- function(lSe)
    {
        lapply(lSe, function(x)
        {
            starts = c(x$starts,x$ends[length(x$ends)])
            s1 = starts[-length(starts)]
            s2 = starts[-c(1)]-1
            nstarts = unlist(lapply(1:length(s1), function(y) round(seq(s1[y],s2[y],3000))))
            nstarts = nstarts[-length(nstarts)]
            nends = c(nstarts[-c(1)]-1,x$ends[length(x$ends)])
            keep = !(nends-nstarts<=0)
            nstarts = nstarts[keep]
            nends = nends[keep]
            list(starts=nstarts,
                 ends=nends)
        })
    }
    lSe <- split_lSe(lSe)
    print("## get all tracks from normal bams")
    require(parallel)
    allTracks <- mclapply(normal_bams,
                          function(bamfile)
                          {
                              cat(".")
                              lCTS.tumour <- lapply(allchr,
                                                    function(chr) getCoverageTrack(bamPath = bamfile,
                                                                                   chr = chr,
                                                                                   lSe[[chr]]$starts,
                                                                                   lSe[[chr]]$ends,
                                                                                   mapqFilter = 30))
                          },
                          mc.cores = MC.CORES)
    print("## calculate normalised ranks")
    rr <- do.call("cbind",lapply(allTracks,function(x) unlist(lapply(x,function(y) y$records))))
    rranks <- apply(rr,2,rank)
    rranks.n <- rranks
    NN <- nrow(rranks)
    for(i in 1:ncol(rranks)) rranks.n[,i] <- rranks.n[,i]/(NN+1-rranks.n[,i])
    rs <- rowSums(log(rranks.n))
    llS <- sapply(1:length(lSe),function(x) length(lSe[[x]]$starts))
    lInds <- lapply(1:length(lSe),function(x) if(x==1) 1:length(lSe[[x]]$starts) else 1:length(lSe[[x]]$starts)+sum(llS[1:(x-1)]))
    ## plot the density of the rank-products and identify the thresholds to remove the tails where outliers are
    plot(density(rs), xlab="ranks across bins", ylab="density", main="")
    title("Choose lower/higher thresholds for outliers")
    ## ---------------------------------##
    ## set thresholds
    ## ---------------------------------##
    cat("Please enter the lower threshold for outliers - a numeric value: ")
    threshold_low <- as.numeric(readline())
    cat("Please enter the higher threshold for outliers - a numeric value: ")
    threshold_high <- as.numeric(readline())
    ## ---------------------------------##
    ## set thresholds
    ## ---------------------------------##
    print("## write and return results")
    lExclude <- lapply(1:length(lSe),function(x)
    {
        cond <- rs[lInds[[x]]]< threshold_low | rs[lInds[[x]]]>threshold_high
        list(starts=lSe[[x]]$starts[cond],
             ends=lSe[[x]]$ends[cond])
    })
    bed <- do.call("rbind",lapply(1:length(lExclude),function(x)
    {
        data.frame(chr=names(lSe)[x],start=lExclude[[x]]$starts,end=lExclude[[x]]$ends)
    }))
    bed$chr=paste0("",bed$chr)
    write.table(bed,file="bed.fasta.targetedexome.bed",
                sep="\t",
                col.names=F,
                row.names=F,
                quote=F)
    return(bed)
}
