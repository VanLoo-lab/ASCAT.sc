printResults_all <- function(res,
                             ismedian=FALSE,
                             outdir="./",
                             is_pdf,
                             projectname="project",
                             svinput=NULL,
                             lSVinput=NULL,
                             rainbowChr=TRUE)
{
    GAMMA <- 1
    if(any(names(res)=="gamma")) GAMMA <- res$gamma
    createDir <- paste0("mkdir ",outdir,"/profiles_",projectname)
    system(createDir)
    for(i in 1:length(res$allTracks.processed))
    {
        if(!is_pdf)
            png(paste0(outdir,"/profiles_",projectname,"/",names(res$allTracks)[i],".png"), width = 5500, height = 2496, res=300)
        if(is_pdf)
            pdf(paste0(outdir,"/profiles_",projectname,"/",names(res$allTracks)[i],".pdf"), width=15, height=5.5)
        try({
            suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                          purity=res$allSolutions[[i]]$purity,
                                          ploidy=res$allSolutions[[i]]$ploidy,
                                          gamma=GAMMA,
                                          ismedian=ismedian,
                                          allchr=res$chr,
                                          ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                          sol=res$allSolutions[[i]],
                                          svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                          rainbowChr=rainbowChr,
                                          is_pdf=is_pdf))
            title(names(res$allTracks)[i])
        })
        dev.off()
        try({
            writeProfile(prof=res$allProfiles[[i]],
                         samplename=paste0(names(res$allTracks)[i],"_",projectname),
                         outdir=outdir)
        })
    }
    ##    zip(zipfile = paste0(outdir,"/profiles_",projectname,".zip"), files = paste0(outdir,"/profiles_",projectname))
    ##    unlink(x=paste0(outdir,"/profiles_",projectname), recursive = TRUE)
    createDirRefit <- paste0("mkdir ",outdir,"/profiles_",projectname,"_refitted")
    system(createDirRefit)
    if(any(grepl("refitted",names(res))))
    {
        for(i in 1:length(res$allTracks.processed))
        {
            if(!is_pdf)
                png(paste0(outdir,"/profiles_",projectname,"_refitted/",names(res$allTracks)[i],".png"), width = 5500, height = 2496, res=300)
            if(is_pdf)
                pdf(paste0(outdir,"/profiles_",projectname,"_refitted/",names(res$allTracks)[i],".pdf"), width = 15, height = 5.5)
            try({
                suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                              purity=res$allSolutions.refitted.auto[[i]]$purity,
                                              ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                              ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                              allchr=res$chr,
                                              gamma=GAMMA,
                                              ismedian=ismedian,
                                              svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                              sol=res$allSolutions[[i]],
                                              rainbowChr=rainbowChr,
                                              is_pdf=is_pdf))
                title(paste0(names(res$allTracks)[i],"-refitted"))
            })
            dev.off()
            try({
                writeProfile(prof=res$allProfiles.refitted.auto[[i]],
                             samplename=paste0(names(res$allTracks)[i],"_",projectname,"_refitted"),
                             outdir=outdir)
            })
        }
        ##zip(zipfile = paste0(outdir,"/profiles_",projectname,"_refitted.zip"), files = paste0(outdir,"/profiles_",projectname,"_refitted"))
        ##unlink(x=paste0(outdir,"/profiles_",projectname,"_refitted"), recursive = TRUE)
    }
    .mytry <- function(x,retVal=NA,...)
    {
        out <- try(x,silent=T,...)
        if(inherits(out,"try-error")) return(retVal)
        out
    }
    getploidy <- function(tt)
    {
        tt <- data.frame(chromosome=as.character(tt[,"chromosome"]),
                         start=as.numeric(tt[,"start"]),
                         end=as.numeric(tt[,"end"]),
                         total_copy_number=as.numeric(tt[,"total_copy_number"]))
        sizes <- (tt$end-tt$start)/1000000
        isna <- is.na(sizes) | is.na(tt$total_copy_number)
        sum(tt$total_copy_number[!isna]*sizes[!isna],na.rm=T)/sum(sizes[!isna],na.rm=T)
    }
    try({res <- append(res,
                       list(summary=list(allSols=data.frame(samplename=names(res$allTracks),
                                                            purity=sapply(res$allSolutions,function(x) .mytry(x$purity)),
                                                            ploidy=sapply(res$allSolutions,function(x) .mytry(x$ploidy)),
                                                            ploidy.tumour=sapply(res$allProfiles,function(x) .mytry(getploidy(x)))),
                                         allSols.refitted=if(!any(grepl("refitted",names(res)))) NULL
                                                          else
                                                              data.frame(samplename=names(res$allTracks),
                                                                         purity=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$purity)),
                                                                         ploidy=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$ploidy)),
                                                                         ploidy.tumour=sapply(res$allProfiles.refitted.auto,function(x) .mytry(getploidy(x)))))))})
    outdir <- gsub("/$","",outdir)
    try(write.table(res$summary$allSols,
                    file=paste0(outdir,
                                "/summary_",
                                projectname,".txt"),
                    sep="\t",quote=F,
                    col.names=T,row.names=T))
    try(write.table(res$summary$allSols.refitted,
                    file=paste0(outdir,"/summary_",projectname,"_refitted.txt"),
                    sep="\t",quote=F,col.names=T,row.names=T))
    save(res, file=paste0(outdir,"/result_object_",projectname,".Rda"))
    res
}
