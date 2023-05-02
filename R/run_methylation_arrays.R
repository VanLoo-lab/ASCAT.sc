run_methylation_array <- function(idat_dir,
                                  res=NULL,
                                  rgSet=NULL,
                                  id_normals=NULL,
                                  sex=NULL,
                                  purs = seq(0.1, 1, 0.01),
                                  ploidies = seq(1.7,5, 0.01),
                                  allchr=c(1:22),
                                  maxtumourpsi=5,
                                  segmentation_alpha=0.0001,
                                  min.width=5,
                                  outdir="./",
                                  conumee=FALSE,
                                  platform=c("450K","Epicv1"),
                                  projectname="project",
                                  predict_refit=TRUE,
                                  print_results=TRUE,
                                  MC.CORES=1)
{
    checkArguments_meth(c(as.list(environment())))
    require(ASCAT.scDataMeth)
    suppressPackageStartupMessages(require(minfi))
    suppressPackageStartupMessages(require(conumee))
    print("## read idat directory")
    .inferSex <- function(tt, annot)
    {
        chrX <- apply(tt[annot[,1]=="X",],2,function(x) median(log2(x),na.rm=T))
        chrY <- apply(tt[annot[,1]=="Y",],2,function(x) median(log2(x),na.rm=T))
        autos <- apply(tt[!annot[,1]%in%c("X","Y"),],2,function(x) median(log2(x),na.rm=T))
        sex <- ifelse(chrX-autos>-.2,ifelse(chrY-autos< -1, "female",NA), ifelse(chrY-autos< -1, NA,"male"))
        sex[is.na(sex)] <- "female" ## default sex
        sex[chrY-autos> -1] <- "male"
        sex
    }
    ## ##################################################
    GAMMA <- .55 ## platform non-linearity parameter - do not change
    anno <- NULL ## initialise bins (not used if conumee=FALSE)
    ## ##################################################
    data("SBDRYs_precomputed",package="ASCAT.sc")
    if(!as.character(segmentation_alpha)%in%names(SBDRYs))
    {
        nperms <- 10000;
        max.ones <- floor(nperms * segmentation_alpha) + 1
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=nperms, max.ones=max.ones)
    }
    else
    {
        SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
    }
    ## ##################################################
    if(is.null(rgSet) & is.null(res))
        suppressWarnings(rgSet <- read.metharray.exp(idat_dir)); gc();
    print("## preprocess raw data - extract unmeth/meth signal unprocessed")
    if(is.null(res))
    {
        data <- preprocessRaw(rgSet); rm("rgSet"); gc();
        ## ##################################################
        totalintensity <- getMeth(data)+getUnmeth(data)+1; gc();
        ## ##################################################
        ## get probe annotation filtered for bad loci
        annot <- as.data.frame(getAnnotation(data)); gc();
        annot[,"chr"] <- gsub("chr","",as.character(annot[,"chr"]))
        annot <- annot[as.character(annot[,"chr"])%in%gsub("chr","",as.character(allchr)),]
        if(platform=="epicv1") data("badloci.epicv1",package="ASCAT.scDataMeth")
        if(platform=="450K")
        {
            data("badloci.450K",package="ASCAT.scDataMeth")
            data("badloci.tcga.450K",package="ASCAT.scDataMeth")
            data("badloci.tcga.450K_segs",package="ASCAT.scDataMeth")
            annot <- annot[!rownames(annot)%in%get(paste0("logr.pon.tcga")) ,]
            annot <- annot[!rownames(annot)%in%get(paste0("bad.loci.tcga_segs")) | annot[,"chr"]%in%c("Y","chrY"),]
        }
        annot <- annot[!rownames(annot)%in%get(paste0("badloci.",platform)) ,]
        ## ##################################################
        ## order chromosome, starts and ends of probes
        annot <- annot[order(annot[,"chr"],annot[,"pos"],decreasing=F),]
        starts <- as.numeric(as.character(annot[,"pos"]))
        ends <- as.numeric(as.character(annot[,"pos"]))
        chrs <- as.character(annot[,"chr"])
        ## ##################################################
        ## Final ordered total intensities for all and normal
        totalintensity <- totalintensity[rownames(annot),]
    }
    if(is.null(id_normals) & !conumee)
    {
        if(is.null(res))
        {
            ## ##################################################
            ## prepare panel of normal (PON) - separate males and females
            if(platform=="450K")
            {
                data("pon.450K",package="ASCAT.scDataMeth")
                data("pon.tcga",package="ASCAT.scDataMeth")
                totalintensityNormal <- cbind(meth_pon[rownames(annot),],pon.tcga[rownames(annot),])
            }
            if(platform=="epicv1")
            {
                data("pon.epicv1",package="ASCAT.scDataMeth")
                totalintensityNormal <- meth_pon[rownames(annot),]
            }
            ## ##################################################
            sexPON <- .inferSex(totalintensityNormal,annot)
            id_normals <- colnames(totalintensityNormal)[sexPON=="female"]
            totalintensityNormalF <- totalintensityNormal[,colnames(totalintensityNormal)%in%id_normals]; gc();
            id_normals <- colnames(totalintensityNormal)[sexPON=="male"]
            totalintensityNormalM <- totalintensityNormal[,colnames(totalintensityNormal)%in%id_normals]; gc();
            ## ##################################################
            if(is.null(sex))
                sex <- .inferSex(totalintensity,annot)
            print("## derive PoN-normalised logr - %%%")
            logr <- sapply(1:ncol(totalintensity),function(x)
            {
                if(sex[x]=="male") totalintensityNormal <- totalintensityNormalM
                if(sex[x]=="female") totalintensityNormal <- totalintensityNormalF
                notalreadyinpanel <- !colnames(totalintensityNormal)%in%colnames(totalintensity)[x]
                TTI <- totalintensity[,x]
                cat(".")
                predicted <- lm(y~.-1,
                                data=data.frame(y=TTI,
                                                X=totalintensityNormal[,notalreadyinpanel]))$fitted.values
                predicted[predicted<1] <- 1
                logr <- log2(totalintensity[,x]/predicted)
                if(sex[x]=="male")
                {
                    logr[annot[,1]=="X"] <- logr[annot[,1]=="X"]-median(logr[annot[,1]%in%c(1:22)])-GAMMA
                    logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-median(logr[annot[,1]%in%c(1:22)])-GAMMA
                }
                if(sex[x]=="female")
                    logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-3.718385
                logr
            })
            cat("\n")
            colnames(logr) <- colnames(totalintensity); gc();
        }
        else
        {
            logr <- res$logr
        }
        print("## segment the non-binned logr tracks")
        allTracks.processed <- mclapply(colnames(logr)[1:ncol(logr)],function(samp)
        {
            cat(".")
            ## ##################################################
            .logr <- meth_winsorise_ascat(logr[,samp])
            ## ##################################################
            getTrackForAll.bins(.logr,
                                paste0("chr",chrs),
                                starts,
                                ends,
                                segmentation_alpha=segmentation_alpha,
                                min.width=min.width,
                                SBDRY=SBDRY,
                                allchr=gsub("chr","",allchr))
            ## ##################################################
        },mc.cores=MC.CORES)
        names(allTracks.processed) <- colnames(logr)[1:ncol(logr)]
        ## ##################################################
    }
    else
    {
        if(length(id_normals)==0 | sum(colnames(totalintensity)%in%id_normals)==0 | is.null(id_normals))
            stop("please provide ids of normal diploid samples or do not use conumee")
        if(is.null(res))
        {
            totalintensityNormal <- totalintensity[,colnames(totalintensity)%in%id_normals]; gc();
            print("## derive PoN-normalised logr")
            ## log PoN-fitted intensity values of all probes for all samples
            logr <- log2(sapply(1:ncol(totalintensity),function(x)
            {
                cat(".")
                notalreadyinpanel <- !colnames(totalintensityNormal)%in%colnames(totalintensity)[x]
                predicted <- lm(y~.-1,
                                data=data.frame(y=totalintensity[,x],
                                                X=totalintensityNormal[,notalreadyinpanel]))$fitted.values
                predicted[predicted<1] <- 1
                totalintensity[,x]/predicted
            }))
            colnames(logr) <- colnames(totalintensity); gc();
            cat("\n")
        }
        else
        {
            logr <- res$logr
        }
        print("## load preprocessed bins dled from https://github.com/mwsill/mnp_training")
        data("CNanalysis4_conumee_ANNO.vh20150715",package="ASCAT.sc")
        print("## bin and segment the binned logr tracks")
        allTracks.processed <- list()
        for(samp in colnames(logr))
        {
            cat(".")
            ## ##################################################
            .logr <- logr[,samp]
            ## ##################################################
            input <- meth_bin(.logr,
                              starts=starts,
                              ends=ends,
                              chrs=paste0("chr",chrs),
                              anno@bins)
            allTracks.processed[[samp]] <- getTrackForAll.bins(input[[1]],
                                                               input[[4]],
                                                               input[[2]],
                                                               input[[3]],
                                                               segmentation_alpha=segmentation_alpha,
                                                               min.width=min.width,
                                                               allchr=gsub("chr","",allchr))
            ## ##################################################
        }
    }
    ## ##################################################
    print("## fit purity and ploidy for all tracks")
    timetofit <- system.time(allSols <- mclapply(1:length(allTracks.processed), function(x)
    {
        sol <- try(searchGrid(allTracks.processed[[x]],
                              purs = purs,
                              ploidies = ploidies,
                              maxTumourPhi=maxtumourpsi,
                              gamma=GAMMA,
                              ismale=sex[x]=="male"),silent=F)
    },mc.cores=MC.CORES))
    ## ##################################################
    print("## get all fitted cna profiles")
    allProfiles <- mclapply(1:length(allTracks.processed), function(x)
    {
        try(getProfile(fitProfile(allTracks.processed[[x]],
                                  purity=allSols[[x]]$purity,
                                  ploidy=allSols[[x]]$ploidy,
                                  gamma=GAMMA,
                                  ismale=sex[x]=="male"),
                       CHRS=allchr),silent=F)
    },mc.cores=MC.CORES)
    names(allProfiles) <- names(allSols) <- names(allTracks.processed)
    ## ##################################################
    print("## refit, print and return all results")
    res <- list(allTracks.processed=allTracks.processed,
                allSolutions=allSols,
                allProfiles=allProfiles,
                chr=allchr,
                sex=sex,
                logr=logr,
                segmentation_alpha=segmentation_alpha,
                min.width=min.width,
                annotions.probes=annot,
                bins=anno,
                gamma=GAMMA,
                timetofit=timetofit)
    if(predict_refit)
        res <- predictRefit_all(res, gamma=GAMMA)
    if(print_results)
        res <- printResults_all(res,  outdir=outdir, projectname=projectname)
    ## ##################################################
    res
}
