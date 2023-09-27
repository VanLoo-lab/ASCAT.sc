printResults_all <- function(res,
                             outdir="./",
                             projectname="project")
{
    GAMMA <- 1
    if(any(names(res)=="gamma")) GAMMA <- res$gamma
    pdf(paste0(outdir,"/profiles_",projectname,".pdf"),width=15,height=4)
    for(i in 1:length(res$allTracks.processed))
    {
        try({
            plotSolution(res$allTracks.processed[[i]],
                         purity=res$allSolutions[[i]]$purity,
                         ploidy=res$allSolutions[[i]]$ploidy,
                         gamma=GAMMA,
                         ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                         sol=res$allSolutions[[i]])
            title(names(res$allTracks)[i])
        })
        try({
            writeProfile(prof=res$allProfiles[[i]],
                         samplename=paste0(names(res$allTracks)[i],"_",projectname),
                         outdir=outdir)
        })
    }
    dev.off()
    if(any(grepl("refitted",names(res))))
    {
        pdf(paste0(outdir,"/profiles_",projectname,"_refitted.pdf"),width=15,height=4)
        for(i in 1:length(res$allTracks.processed))
        {
            try({
                plotSolution(res$allTracks.processed[[i]],
                             purity=res$allSolutions.refitted.auto[[i]]$purity,
                             ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                             ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                             gamma=GAMMA,
                             sol=res$allSolutions[[i]])
                title(paste0(names(res$allTracks)[i],"-refitted"))
            })
            try({
                writeProfile(prof=res$allProfiles.refitted.auto[[i]],
                             samplename=paste0(names(res$allTracks)[i],"_",projectname),
                             outdir=outdir)
            })
        }
        dev.off()
    }
    .mytry <- function(x,retVal=NA,...)
    {
        out <- try(x,silent=T,...)
        if(inherits(out,"try-error")) return(retVal)
        out
    }
    try({res <- append(res,
                       list(summary=list(allSols=data.frame(samplename=names(res$allTracks),
                                                            purity=sapply(res$allSolutions,function(x) .mytry(x$purity)),
                                                            ploidy=sapply(res$allSolutions,function(x) .mytry(x$ploidy))),
                                         allSols.refitted=if(!any(grepl("refitted",names(res)))) NULL
                                                          else
                                                              data.frame(samplename=names(res$allTracks),
                                                                         purity=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$purity)),
                                                                         ploidy=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$ploidy))))))})
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
    res
}
