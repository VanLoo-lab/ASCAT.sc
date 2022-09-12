run_methylation_array <- function(idat_dir,
                                  id_normals,
                                  purs = seq(0.1, 1, 0.01),
                                  ploidies = seq(1.7,5, 0.01),
                                  allchr=c(1:22),
                                  maxtumourpsi=5,
                                  segmentation_alpha=0.01,
                                  outdir="./",
                                  projectname="project",
                                  predict_refit=TRUE,
                                  print_results=TRUE,
                                  MC.CORES=1)
{
    suppressPackageStartupMessages(require(minfi))
    suppressPackageStartupMessages(require(conumee))
    print("## read idat directory")
    rgSet <- read.metharray.exp(idat_dir); gc();
    print("## preprocess raw data - extract unmeth/meth signal unprocessed")
    data <- preprocessRaw(rgSet); rm("rgSet"); gc();
    print("## load preprocessed bins dled from https://github.com/mwsill/mnp_training")
    data("CNanalysis4_conumee_ANNO.vh20150715",package="ASCAT.sc")
    ## ##################################################
    totalintensity <- getMeth(data)+getUnmeth(data)+1; gc();
    annot <- as.data.frame(getAnnotation(data)); rm("data"); gc();
    annot <- annot[as.character(annot[,"chr"])%in%as.character(allchr),]
    ## ##################################################
    ## order chromosome, starts and ends of probes
    annot <- annot[order(annot[,"chr"],annot[,"pos"],decreasing=F),]
    starts <- as.numeric(as.character(annot[,"pos"]))
    ends <- as.numeric(as.character(annot[,"pos"]))
    chrs <- as.character(annot[,"chr"])
    ## ##################################################
    ## Final ordered total intensities for all and normal
    totalintensity <- totalintensity[rownames(annot),]
    totalintensityNormal <- totalintensity[,colnames(totalintensity)%in%id_normals]; gc();
    ## ##################################################
    ## ##################################################
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
    ## ##################################################
    cat("\n")
    ## ##################################################
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
                          chrs=chrs,
                          anno@bins)
        allTracks.processed[[samp]] <- getTrackForAll.bins(input[[1]],
                                                           input[[4]],
                                                           input[[2]],
                                                           input[[3]],
                                                           allchr=gsub("chr","",allchr))
        ## ##################################################
    }
    ## ##################################################
    GAMMA <- .55 ## platform non-linearity parameter - do not change
    ## ##################################################
    print("## fit purity and ploidy for all tracks")
    timetofit <- system.time(allSols <- mclapply(1:length(allTracks.processed), function(x)
    {
        sol <- try(searchGrid(allTracks.processed[[x]],
                              purs = purs,
                              ploidies = ploidies,
                              maxTumourPhi=maxtumourpsi,
                              gamma=GAMMA,
                              ismale=F),silent=F)
    },mc.cores=MC.CORES))
    ## ##################################################
    print("## get all fitted cna profiles")
    allProfiles <- mclapply(1:length(allTracks.processed), function(x)
    {
        try(getProfile(fitProfile(allTracks.processed[[x]],
                                  purity=allSols[[x]]$purity,
                                  ploidy=allSols[[x]]$ploidy,
                                  gamma=GAMMA,
                                  ismale=F),
                       CHRS=allchr),silent=F)
    },mc.cores=MC.CORES)
    names(allProfiles) <- names(allSols) <- names(allTracks.processed)
    ## ##################################################
    print("## refit, print and return all results")
    res <- list(allTracks.processed=allTracks.processed,
                allSolutions=allSols,
                allProfiles=allProfiles,
                chr=allchr,
                sex=rep("female",length(allTracks.processed)),
                logr=logr,
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
