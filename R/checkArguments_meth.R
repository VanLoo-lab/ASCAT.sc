checkArguments_meth <- function(args)
{
    print("checking arguments")
    nms <- c("allTracks.processed",
             "allSolutions",
             "allProfiles",
             "chr","sex","logr","segmentation_alpha",
             "min.width","annotations.probes","bins","gamma")
    nfiles <- length(dir(args$idat_dir,pattern="idat$"))
    if(nfiles==0 & !is.null(args$idat_dir)) stop("No files found in directory")
    if(!is.null(args$res))
    {
        if(!all(nms%in%names(args$res)))
<<<<<<< HEAD
            stop(paste0("Malformed res object - missing elements ",nms[!nms%in%names(args$res)]))
=======
            stop(paste0("Malformed res object - missing elements",nms[!nms%in%names(args$res)]))
>>>>>>> d619bd95c4c64791b532ea98ba54d68877ef031b
    }
    if(!dir.exists(args$outdir))
        stop("Output directory outdir does not exist")
}
