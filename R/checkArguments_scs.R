checkArguments_scs <- function(args)
{
    print("checking arguments")
    nms <- c("allTracks.processed",
             "allTracks",
             "allSolutions",
             "allProfiles",
             "purs","ploidies","maxtumourpsi",
             "chr","sex","logr","segmentation_alpha",
             "min.width","annotations.probes","bins","gamma")
    if(!is.null(args$res))
        if(!all(nms%in%names(args$res)))
            stop(paste0("Malformed res object - missing elements",nms[!nms%in%names(args$res)]))
    if(!dir.exists(args$outdir))
        stop("Output directory outdir does not exist")
    if(length(args$sex)!=length(args$tumour_bams))
        print("Note: not provided a sex (female or male) for each bam")
}
