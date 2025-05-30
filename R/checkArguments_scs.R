checkArguments_scs <- function(args)
{
    print("checking arguments")
    nms <- c("allTracks.processed",
             "allTracks",
             "allSolutions",
             "allProfiles",
             "purs",
             "ploidies","maxtumourpsi",
             "chr","sex",
             "segmentation_alpha")
    if(!is.null(args$res))
        if(!all(nms%in%names(args$res)))
            stop(paste("Malformed res object - missing elements",nms[!nms%in%names(args$res)],collapse=" "))
    if(!dir.exists(args$outdir))
        stop("Output directory outdir does not exist")
    if(length(args$sex)!=length(args$tumour_bams))
        print("Note: not provided a sex (female or male) for each bam")
}
