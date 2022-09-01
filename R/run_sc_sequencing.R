run_sc_sequencing <- function(tumour_bams,
                              res=NULL,
                              allchr=paste0("",c(1:22,"X")),
                              sex=c("female","male","female"),
                              chrstring_bam="",
                              purs = seq(0.1, 1, 0.01),
                              ploidies = seq(1.7,5, 0.01),
                              maxtumourpsi=5,
                              binsize=500000,
                              segmentation_alpha=0.01,
                              predict_refit=TRUE,
                              print_results=TRUE,
                              build=c("hg19","hg38"),
                              MC.CORES=1,
                              barcodes_10x=NULL,
                              outdir="./",
                              probs_filters=.1,
                              path_to_phases=NULL,
                              list_ac_counts_paths=NULL,
                              sc_filters=FALSE,
                              projectname="project",
                              smooth_sc=FALSE,
                              multipcf=FALSE)
{
    require(parallel)
    require(Rsamtools)
    require(Biostrings)
    require(DNAcopy)
    require(copynumber)
    binsize <- as.numeric(binsize)
    print(paste0("## load Bins for Genome Build: ", build))
    START_WINDOW <- 30000
    if(is.null(res))
        res <- list()
    if(build=="hg19")
    {
        data("lSe_filtered_30000.hg19",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg19",package="ASCAT.sc")
        if(!all(allchr%in%names(lSe.hg19.filtered)))
            stop(paste0("allchr should be in the form: ",names(lSe.hg19.filtered)[1],"; not: ",allchr[1],"..."))
        res$lSe <- lapply(allchr, function(chr) lSe.hg19.filtered[[chr]])
        names(res$lSe) <- allchr
        names(lGCT.hg19.filtered) <- names(res$lSe)
        res$lGCT <- lapply(allchr, function(chr) lGCT.hg19.filtered[[chr]])
    }
    if(build=="hg38")
    {
        data("lSe_filtered_30000.hg38",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg38",package="ASCAT.sc")
        res$lSe <- lapply(allchr, function(chr) lSe.hg38.filtered[[chr]])
        names(res$lSe) <- allchr
        names(lGCT.hg38.filtered) <- names(res$lSe)
        res$lGCT <- lapply(allchr, function(chr) lGCT.hg38.filtered[[chr]])
        names(res$lGCT) <- names(res$lSe)
        if(chrstring_bam=="")
            names(res$lGCT) <- names(res$lSe) <- gsub("chr","",names(res$lSe))
    }
    print("## calculate Target Bin size")
    res$nlGCT <- treatGCT(res$lGCT,window=ceiling(binsize/START_WINDOW))
    res$nlSe <- treatlSe(res$lSe,window=ceiling(binsize/START_WINDOW))
    if(is.null(res) | !any(names(res)=="allTracks"))
    {
        if(!is.null(barcodes_10x))
        {
            print("## get Tracks from 10X-like bam")
            timetoread_tumours <- system.time(res$allTracks <- getTrackForAll.10XBAM(bamfile=tumour_bams[1],
                                                                                     barcodes=barcodes_10x,
                                                                                     allchr=allchr,
                                                                                     chrstring="",
                                                                                     lSe=res$lSe,
                                                                                     doSeg=FALSE))
            nms <- names(res$allTracks)
            res$allTracks <- lapply(nms, function(x)
                list(lCTS.tumour=res$allTracks[[x]],
                     nlCTS.tumour=treatTrack(lCTS=res$allTracks[[x]],
                                             window=ceiling(binsize/START_WINDOW)))
                )
            names(res$allTracks) <- nms
            sex <- rep(sex,length(res$allTracks))
        }
        if(is.null(barcodes_10x))
        {
            print("## get Tracks from Tumour bams")
            timetoread_tumours <- system.time(res$allTracks <- mclapply(tumour_bams,function(bamfile)
            {
                lCTS.tumour <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                                             chr=chr,
                                                                             res$lSe[[chr]]$starts,
                                                                             res$lSe[[chr]]$ends,
                                                                             mapqFilter=30))
                list(lCTS.tumour=lCTS.tumour,
                     nlCTS.tumour=treatTrack(lCTS=lCTS.tumour,
                                             window=ceiling(binsize/START_WINDOW)))
            },mc.cores=MC.CORES))
        }
    }
    if(multipcf)
    {
        print("## calculate Multipcf - multi-sample mode - do not use if samples from different tumours")
        timetoprocessed <- system.time(res$allTracks.processed <- getLSegs.multipcf(allTracks=lapply(res$allTracks, function(x)
                                                                                                {list(lCTS=x$nlCTS.tumour)}),
                                                                                lCTS=lapply(res$allTracks,function(x) x$nlCTS.tumour),
                                                                                lSe=res$nlSe,
                                                                                lGCT=res$nlGCT,
                                                                                lNormals=NULL,
                                                                                allchr=allchr,
                                                                                segmentation_alpha=segmentation_alpha,
                                                                                MC.CORES=MC.CORES))
    }
    else
    {
        print("## smooth Tracks")
        timetoprocessed <- system.time(res$allTracks.processed <- mclapply(1:length(res$allTracks), function(x)
        {
            getTrackForAll(bamfile=NULL,
                           window=NULL,
                           lCT=res$allTracks[[x]][[2]],
                           lSe=res$nlSe,
                           lGCT=res$nlGCT,
                           lNormals=NULL,
                           allchr=allchr,
                           sdNormalise=0,
                           segmentation_alpha=segmentation_alpha)
        },mc.cores=MC.CORES))
        cat("\n")
    }
    if(is.null(barcodes_10x))
        names(res$allTracks.processed) <- names(res$allTracks) <- gsub("(.*)/(.*)","\\2",tumour_bams)
    else
        names(res$allTracks.processed) <- names(res$allTracks)
    print("## fit Purity/Ploidy")
    timetofit <- system.time(res$allSols <- mclapply(1:length(res$allTracks.processed), function(x)
    {
        sol <- try(searchGrid(res$allTracks.processed[[x]],
                              purs = purs,
                              ploidies = ploidies,
                              maxTumourPhi=maxtumourpsi,
                              ismale=if(sex[x]=="male") T else F),silent=F)
    },mc.cores=MC.CORES))
    print("## get Fitted Cna Profiles")
    res$allProfiles <- mclapply(1:length(res$allTracks.processed), function(x)
    {
        try(getProfile(fitProfile(res$allTracks.processed[[x]],
                                  purity=res$allSols[[x]]$purity,
                                  ploidy=res$allSols[[x]]$ploidy,
                                  ismale=if(sex[x]=="male") T else F),
                       CHRS=allchr),silent=F)
    },mc.cores=MC.CORES)
    names(res$allProfiles) <- names(res$allSols) <- names(res$allTracks)
    print("## compile Results")
    res <- list(allTracks.processed=res$allTracks.processed,
                allTracks=res$allTracks,
                allSolutions=res$allSols,
                allProfiles=res$allProfiles,
                chr=paste0(chrstring_bam,allchr),
                purs=purs,
                ploidies=ploidies,
                maxtumourpsi=maxtumourpsi,
                build=build,
                binsize=binsize,
                sex=sex,
                segmentation_alpha=segmentation_alpha,
                multipcf=multipcf,
                lSe=res$nlSe,
                lGCT=res$nlGCT,
                timetoread_tumours=timetoread_tumours,
                timetoprocessed=timetoprocessed,
                timetofit=timetofit)
    if(predict_refit)
    {
        print("## predict Refit")
        res <- predictRefit_all(res)
    }
    if(print_results)
    {
        print("## print Results")
        res <- printResults_all(res, outdir=outdir, projectname=projectname)
    }
    if(sc_filters)
    {
        print("## get Filters")
        res<- getFilters(res,
                         probs=probs_filters,
                         outdir=outdir,
                         projectname=projectname)
    }
    if(!is.null(list_ac_counts_paths) & !is.null(path_to_phases))
    {
        print("## get Allele-specific CNA")
        res <- getAS_CNA(res,
                         path_to_phases=path_to_phases,
                         list_ac_counts_paths=list_ac_counts_paths,
                         purs=purs,
                         ploidies=ploidies,
                         outdir=outdir,
                         projectname=projectname,
                         mc.cores=MC.CORES)
        names(res$allProfiles_AS) <- names(res$allProfiles)
    }
    if(smooth_sc & any(grepl("_AS",names(res))))
    {
        print("## smooth Across Single-Cells/Nulcei")
        res <- getAS_CNA_smoothed(res,
                                  mc.cores=MC.CORES)
    }
    if(smooth_sc & !any(grepl("_AS",names(res))))
        print("Warning: Smoothing is only possible for allele-specific copy numbers")
    res
}
