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
                              pcchromosome = 0.8,
                              predict_refit=TRUE,
                              print_results=TRUE,
                              build=c("hg19","hg38"),
                              MC.CORES=1,
                              barcodes_10x=NULL,
                              normal_bams=NULL,
                              outdir="./",
                              probs_filters=.1,
                              path_to_phases=NULL,
                              list_ac_counts_paths=NULL,
                              sc_filters=FALSE,
                              projectname="project",
                              steps=NULL,
                              smooth_sc=FALSE,
                              multipcf=FALSE,
                              normalize=FALSE,
                              lExclude=NULL,
                              sc_exclude_badbins=FALSE)
{
    checkArguments_scs(c(as.list(environment())))
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(Biostrings))
    suppressPackageStartupMessages(require(DNAcopy))
    suppressPackageStartupMessages(require(copynumber))
    binsize <- as.numeric(binsize)
    print(paste0("## load Bins for Genome Build: ", build))
    START_WINDOW <- 30000
    data("SBDRYs_precomputed",package="ASCAT.sc")
    if(!as.character(segmentation_alpha)%in%names(SBDRYs))
    {
        nperms <- 10000;
        max.ones <- floor(nperms * segmentation_alpha) + 1
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=nperms, max.ones=max.ones)
    } else
    {
        SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
    }
    if(binsize<30000)
    {
        print("Current minimum bin size is 30000 - resetting to 30000")
    }
    if(is.null(res))
        res <- list()
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
        allchr. <- paste0("chr",gsub(chrstring_bam,"",allchr))
        res$lSe <- lapply(allchr., function(chr) lSe.hg38.filtered[[chr]])
        names(res$lSe) <- allchr
        names(lGCT.hg38.filtered) <- names(res$lSe)
        res$lGCT <- lapply(allchr., function(chr) lGCT.hg38.filtered[[chr]])
        names(res$lGCT) <- names(res$lSe)
        if(chrstring_bam=="")
            names(res$lGCT) <- names(res$lSe) <- gsub("chr","",names(res$lSe))
    }
    print("## calculate Target Bin size")
    res$nlGCT <- treatGCT(res$lGCT,window=ceiling(binsize/START_WINDOW))
    res$nlSe <- treatlSe(res$lSe,window=ceiling(binsize/START_WINDOW))
    if(is.null(res$isPON)) res$isPON <- FALSE
    if(!any(names(res)=="binsize")) res$binsize <- binsize
    if(!is.null(barcodes_10x))
    {
        if(!any(names(res)=="allTracks"))
        {
            print("## get Tracks from 10X-like bam")
            res$timetoread_tumours <- system.time(res$allTracks <- lapply(tumour_bams,function(bamfile)
            {
                lCTS.tumour <- getTrackForAll.10XBAM(bamfile=bamfile,
                                                     barcodes=barcodes_10x,
                                                     allchr=allchr,
                                                     chrstring="",
                                                     lSe=res$lSe,
                                                     doSeg=FALSE,
                                                     # Percent of chromosome bins that need to have a barcode to be kept
                                                     pcchromosome = pcchromosome, 
                                                     mc.cores=MC.CORES)

                # If multiple bams, rename the tracks to include bamfile (sample)
                
                if (length(tumour_bams) > 1) {
                    names(lCTS.tumour) <- paste0(basename(bamfile), "_", names(lCTS.tumour))
                }

                # Fix binsizes to match the desired binsize using treatTrack()
                mclapply(lCTS.tumour, function(x) {
                    list(lCTS.tumour = x,
                         nlCTS.tumour = treatTrack(lCTS = x,
                                                   window = ceiling(binsize / START_WINDOW)))
                }, mc.cores = MC.CORES)
            }))
            # Unlist the top layer to flatten list
            res$allTracks <- unlist(res$allTracks, recursive = FALSE)

            # sex[1] because the assumption is when running multiple bamfiles, they are from the same tumour?
            # This ensures that it has the correct length
            sex <- rep(sex[1], length(res$allTracks))
        }
    }
    if(is.null(barcodes_10x))
    {
        if(!is.null(normal_bams[1]) & is.null(res$nlCTS.normal))
        {
            print("## get all tracks from normal bams")
            timetoread_normals <- system.time(res$lCTS.normal <- mclapply(normal_bams,function(bamfile)
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
            res$nlCTS.normal <- combineDiploid(lapply(res$lCTS.normal,function(x) x[[2]]))
            res$lNormals <- lapply(res$lCTS.normal,function(x) x$nlCTS.normal)
            res$isPON <- TRUE
        }
        else
        {
            res$isPON <- FALSE
            res$nlCTS.normal <- NULL
            res$lNormals <- NULL
            res$timetoread_normals <- NULL
        }
        if(any(names(res)=="allTracks") & res$binsize!=binsize)
        {
            print("## adjust Tracks for bin size ")
            res$timetoread_tumours <- system.time(res$allTracks <- mclapply(names(res$allTracks),function(bamfile)
            {
                lCTS.tumour <- res$allTracks[[bamfile]]$lCTS.tumour
                list(lCTS.tumour=lCTS.tumour,
                     nlCTS.tumour=treatTrack(lCTS=lCTS.tumour,
                                             window=ceiling(binsize/START_WINDOW)))
            },mc.cores=MC.CORES))
        }
        if(!any(names(res)=="allTracks"))
        {
            print("## get Tracks from Tumour bams")
            res$timetoread_tumours <- system.time(res$allTracks <- mclapply(tumour_bams,function(bamfile)
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
            names(res$allTracks) <- basename(tumour_bams)
        }
    }
    if(sc_exclude_badbins)
    {
        print("## exclude Bad bins inferred from normal samples")
        res <- sc_excludeBadBins(res)
    }
    if(multipcf)
    {
        print("## calculate Multipcf - multi-sample mode - do not use if samples from different tumours")
        res$timetoprocessed <- system.time(res$allTracks.processed <- getLSegs.multipcf(allTracks=lapply(res$allTracks, function(x)
        {list(lCTS=x$nlCTS.tumour)}),
        lCTS=lapply(res$allTracks,function(x) x$nlCTS.tumour),
        lSe=res$nlSe,
        lGCT=res$nlGCT,
        lNormals=res$lNormals,
        allchr=allchr,
        segmentation_alpha=segmentation_alpha,
        normalize=normalize,
        MC.CORES=MC.CORES))
    }
    else
    {
        if(!any("allTracks.processed"%in%names(res)) | res$binsize!=binsize)
        {
            print("## smooth Tracks")
            res$timetoprocessed <- system.time(res$allTracks.processed <- mclapply(1:length(res$allTracks), function(x)
            {
                getTrackForAll(bamfile=NULL,
                               window=NULL,
                               lCT=res$allTracks[[x]][[2]],
                               lSe=res$nlSe,
                               lGCT=res$nlGCT,
                               lNormals=res$lNormals,
                               allchr=allchr,
                               sdNormalise=0,
                               SBDRY=SBDRY,
                               segmentation_alpha=segmentation_alpha)
            },mc.cores=MC.CORES))
        }
    }
    if(is.null(barcodes_10x))
        names(res$allTracks.processed) <- names(res$allTracks) <- basename(tumour_bams)
    else
        names(res$allTracks.processed) <- names(res$allTracks)
    print("## fit Purity/Ploidy")
    res$timetofit <- system.time(res$allSols <- mclapply(1:length(res$allTracks.processed), function(x)
    {
        sol <- try(searchGrid(res$allTracks.processed[[x]],
                              purs = purs,
                              ploidies = ploidies,
                              maxTumourPhi=maxtumourpsi,
                              ismale=if(sex[x]=="male") T else F,
                              isPON=res$isPON),silent=F)
    },mc.cores=MC.CORES))
    print("## get Fitted CNA Profiles")
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
                chr=allchr,
                chrstring_bam=chrstring_bam,
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
                isPON=res$isPON,
                timetoread_tumours=res$timetoread_tumours,
                timetoprocessed=res$timetoprocessed,
                timetofit=res$timetofit,
                mode="sc")
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
                         steps=steps,
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
