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
            stop(paste0("Malformed res object - missing elements ",nms[!nms%in%names(args$res)]))
    }
    if(!dir.exists(args$outdir))
        stop("Output directory outdir does not exist")
    getPlatform <- function(file)
    {
        require(minfi)
        tryCatch({
            readIDAT(file, quiet = TRUE)$RunInfo$Array
        }, error = function(e) NA)
    }
    idat_files <- list.files(args$idat_dir, pattern = "idat$", full.names = TRUE)
    platforms <- sapply(idat_files, getPlatform)
    print(table(platforms))

}
