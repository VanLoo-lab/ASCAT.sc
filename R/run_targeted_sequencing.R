run_targeted_sequencing <- function(tumour_bams,
                                    bed_file,
                                    allchr=paste0("",c(1:22,"X")),
                                    sex=c("female","male","female"),
                                    chrstring_bed="",
                                    chrstring_bam="",
                                    purs = seq(0.1, 1, 0.01),
                                    ploidies = seq(1.7,5, 0.01),
                                    maxtumourpsi=5,
                                    bed_padding=1000,
                                    binsize=500000,
                                    segmentation_alpha=0.01,
                                    normal_bams=NULL,
                                    nlCTS.normal=NULL,
                                    build=c("hg19","hg38"),
                                    predict_refit=TRUE,
                                    print_results=TRUE,
                                    MC.CORES=1,
                                    outdir="./",
                                    projectname="project",
                                    multipcf=FALSE)
{
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(Biostrings))
    suppressPackageStartupMessages(require(DNAcopy))
    suppressPackageStartupMessages(require(copynumber))
    binsize <- as.numeric(binsize)
    print("## load bins for genome build")
    START_WINDOW <- 30000
    if(build=="hg19")
    {
        data("lSe_filtered_30000.hg19",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg19",package="ASCAT.sc")
        if(!all(allchr%in%names(lSe.hg19.filtered)))
            stop(paste0("allchr should be in the form: ",names(lSe.hg19.filtered)[1],"; not: ",allchr[1],"..."))
        lSe <- lapply(allchr, function(chr) lSe.hg19.filtered[[chr]])
        names(lGCT.hg19.filtered) <- names(lSe.hg19.filtered)
        names(lSe) <- paste0(chrstring_bam,allchr)
        lGCT <- lapply(allchr, function(chr) lGCT.hg19.filtered[[chr]])
        names(lGCT) <- names(lSe)
    }
    if(build=="hg38")
    {
        data("lSe_filtered_30000.hg38",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg38",package="ASCAT.sc")
        lSe <- lapply(allchr, function(chr) lSe.hg38.filtered[[chr]])
        names(lSe) <- allchr
        names(lGCT.hg38.filtered) <- names(lSe)
        lGCT <- lapply(allchr, function(chr) lGCT.hg38.filtered[[chr]])
        names(lGCT) <- names(lSe)
        if(chrstring_bam=="")
            names(res$lGCT) <- names(res$lSe) <- gsub("chr","",names(res$lSe))
    }
    print("## read in bed file")
    bed <- ts_treatBed(read.table(bed_file,
                                  sep="\t",
                                  header=F),
                       add=bed_padding)
    print("## get bin starts ends to exclude")
    lExclude <- lapply(allchr,function(chr)
        ts_getExcludeFromBedfile(bed,chr))
    names(lExclude) <- allchr
    lInds <- getlInds(lSe, lExclude)
    lGCT <- getlGCT_excluded(lGCT, lInds)
    lSe <- getlSe_excluded(lSe, lInds)
    lCTS.normal <- NULL
    lCTS.normal.combined <- NULL
    timetoread_normals <- NULL
    if(!is.null(normal_bams[1]) & is.null(nlCTS.normal))
    {
        print("## get all tracks from normal bams")
        timetoread_normals <- system.time(lCTS.normal <- mclapply(normal_bams,function(bamfile)
        {
            lCTS.normal <- lapply(paste0(chrstring_bam,allchr), function(chr) getCoverageTrack(bamPath=bamfile,
                                                                         chr=chr,
                                                                         lSe[[chr]]$starts,
                                                                         lSe[[chr]]$ends,
                                                                         mapqFilter=30))
            list(lCTS.normal=lCTS.normal,
                 nlCTS.normal=treatTrack(lCTS=lCTS.normal,
                                        window=ceiling(binsize/START_WINDOW)))
        },mc.cores=MC.CORES))
        lCTS.normal.combined <- combineDiploid(lapply(lCTS.normal,function(x) x[[2]]))
    }
    print("## calculate target bin size")
    nlGCT <- treatGCT(lGCT,window=ceiling(binsize/START_WINDOW))
    nlSe <- treatlSe(lSe,window=ceiling(binsize/START_WINDOW))
    print("## get all tracks from tumour bams")
    timetoread_tumours <- system.time(allTracks <- mclapply(tumour_bams,function(bamfile)
    {
        lCTS.tumour <- lapply(paste0(chrstring_bam,allchr), function(chr) getCoverageTrack(bamPath=bamfile,
                                                                     chr=chr,
                                                                     lSe[[chr]]$starts,
                                                                     lSe[[chr]]$ends,
                                                                     mapqFilter=30))
        list(lCTS.tumour=lCTS.tumour,
             nlCTS.tumour=treatTrack(lCTS=lCTS.tumour,
                                     window=ceiling(binsize/START_WINDOW)))
    },mc.cores=MC.CORES))
    if(multipcf)
    {
        print("## calculating multipcf - multi-sample mode - do not use if samples from different tumours")
        timetoprocessed <- system.time(allTracks.processed <- getLSegs.multipcf(allTracks=lapply(allTracks, function(x)
                                                                                                {list(lCTS=x$nlCTS.tumour)}),
                                                                                lCTS=lapply(allTracks,function(x) x$nlCTS.tumour),
                                                                                lSe=nlSe,
                                                                                lGCT=nlGCT,
                                                                                lNormals=lNormals,
                                                                                allchr=allchr,
                                                                                segmentation_alpha=segmentation_alpha,
                                                                                MC.CORES=MC.CORES))
    }
    else
    {
        print("## smooth and (apply) segments to all tracks")
        timetoprocessed <- system.time(allTracks.processed <- mclapply(1:length(allTracks), function(x)
        {
            cat(".")
            getTrackForAll(bamfile=NULL,
                           window=NULL,
                           lCT=allTracks[[x]][[2]],
                           lSe=nlSe,
                           lGCT=nlGCT,
                           lNormals=lCTS.normal.combined,
                           allchr=allchr,
                           sdNormalise=0,
                           segmentation_alpha=segmentation_alpha)
        },mc.cores=MC.CORES))
        cat("\n")
    }
    names(allTracks.processed) <- names(allTracks) <- gsub("(.*)/(.*)","\\2",tumour_bams)
    print("## fit purity and ploidy for all tracks")
    timetofit <- system.time(allSols <- mclapply(1:length(allTracks.processed), function(x)
    {
        sol <- try(searchGrid(allTracks.processed[[x]],
                              purs = purs,
                              ploidies = ploidies,
                              maxTumourPhi=maxtumourpsi,
                              ismale=if(sex[x]=="male") T else F),silent=F)
    },mc.cores=MC.CORES))
    print("## get all fitted cna profiles")
    allProfiles <- mclapply(1:length(allTracks.processed), function(x)
    {
        try(getProfile(fitProfile(allTracks.processed[[x]],
                                  purity=allSols[[x]]$purity,
                                  ploidy=allSols[[x]]$ploidy,
                                  ismale=if(sex[x]=="male") T else F),
                       CHRS=paste0(chrstring_bam,allchr)),silent=F)
    },mc.cores=MC.CORES)
    names(allProfiles) <- names(allSols) <- names(allTracks)
    print("## return all results")
    res <- list(allTracks.processed=allTracks.processed,
                allTracks=allTracks,
                allSolutions=allSols,
                allProfiles=allProfiles,
                lCTS.normal=lCTS.normal,
                nlCTS.normal=nlCTS.normal,
                chr=allchr,
                chrstring_bam=chrstring_bam,
                sex=sex,
                lSe=nlSe,
                purs=purs,
                ploidies=ploidies,
                maxtumourpsi=maxtumourpsi,
                lExclude=lExclude,
                multipcf=multipcf,
                lGCT=nlGCT,
                build=build,
                binsize=binsize,
                timetoread_normals=timetoread_normals,
                timetoread_tumours=timetoread_tumours,
                timetoprocessed=timetoprocessed,
                timetofit=timetofit)
    if(predict_refit)
        res <- predictRefit_all(res)
    if(print_results)
        res <- printResults_all(res, outdir=outdir, projectname=projectname)
    res
}
