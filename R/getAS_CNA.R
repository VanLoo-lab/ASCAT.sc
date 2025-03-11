getAS_CNA <- function(res,
                      path_to_phases,
                      list_ac_counts_paths,
                      purs,
                      ploidies,
                      chrstring="chr",
                      projectname="project",
                      outdir="./",
                      mc.cores=1,
                      steps=NULL)
{
    suppressPackageStartupMessages(require(GenomicRanges))
    suppressPackageStartupMessages(require(data.table))
    getNANBfromTot <- function(baf,tot,purity, retNa=TRUE)
    {
        Na <- round(baf*(purity*tot+(1-purity)*2)-1+purity)
        Nb <- round(tot-Na)
        if(retNa) return(Na)
        return(Nb)
    }
    .searchGrid  <-  function (baf,
                              logr,
                              sizes,
                              purs=seq(.90,.99,.01),
                              ploidies=seq(2,6,.01))
    {
        getNANB <- function(baf,logr,purity,ploidy)
        {
            K <- (2^logr*ploidy-2*(1-purity))/purity
            Na <- (1-purity-baf*(2-2*purity)-(purity*baf-purity)*K)/purity
            Nb <- K-Na
            list(Na=Na,Nb=Nb)
        }

        errors <- function(baf, logr, sizes, purity, ploidy)
        {
            nanb <- getNANB(baf,logr,purity,ploidy)
            ssize <- sum(sizes)
            sum(((nanb$Na-round(nanb$Na))^2+(nanb$Nb-round(nanb$Nb))^2)*(sizes/ssize))
        }

        ggprofile <- function(baf, logr, purity, ploidy)
        {
            nanb <- getNANB(baf,logr,purity,ploidy)
        }
        errs <- t(sapply(purs,function(rho)
            sapply(ploidies,function(psi)
                errors(baf=baf,
                       logr=logr,
                       sizes,
                       purity=rho,ploidy=psi))))
        rownames(errs) <- purs
        colnames(errs) <- ploidies
        purs <- as.numeric(rownames(errs))
        ploidies <- as.numeric(colnames(errs))
        mins <- arrayInd(which.min(errs), dim(errs))
        purity <- purs[mins[1]]
        ploidy <- ploidies[mins[2]]
        return(list(purity=purity,
                    ploidy=ploidy,
                    profile=ggprofile(baf,logr,purity,ploidy),
                    errs=errs))
    }

    readPhases <- function(phasing_paths)
    {
        phasing <- lapply(phasing_paths,function(x) as.data.frame(data.table::fread(x)))
        phases <- lapply(phasing,function(x)
        {
            ##x <- x[x[,10]%in%c("0|1","1|0"),]
            x <- x[grep("0\\|1|1\\|0",x[, 10]), ]
            phase <- gsub("(.*)\\|(.*)","\\1",x[,10])
            phases1 <- x[,"REF"]
            phases1[phase=="1"] <- x[phase=="1","ALT"]
            phases2 <- x[,"ALT"]
            phases2[phases2==phases1] <- x[phases2==phases1,"REF"]
            list(chr=x[,1],
                 pos=x[,2],
                 phases1=phases1,
                 phases2=phases2)
        })
        phasesall <- do.call("rbind",lapply(phases,function(x) data.frame(chr=x[[1]],
                                                                          pos=x[[2]],
                                                                          phase1=x[[3]],
                                                                          phase2=x[[4]])))
    }

    getPhasedInfo <- function(ac, phasing)
    {
        ac[,1] <- gsub("chr","",ac[,1])
        phasing[,1] <- gsub("chr","",phasing[,1])
        ac.. <- merge(ac, phasing,by.x=c("#CHR","POS"), by.y=c("chr","pos"))
        nmsPH1 <- ac..[,"phase1"]
        nmsPH2 <- ac..[,"phase2"]
        letters <- c("A","C","G","T")
        letters_index <- sapply(paste0("Count_",letters),function(x) which(colnames(ac..)%in%x))
        names(letters_index) <-letters
        ind1 <- letters_index[nmsPH1]
        ind2 <- letters_index[nmsPH2]
        ind1 <- cbind(seq_along(ind1), ind1)
        ind2 <- cbind(seq_along(ind2), ind2)
        df <- data.frame(chr=ac..[,1],
                         pos=ac..[,2],
                         counts1=as.numeric(as.character(ac..[ind1])),
                         counts2=as.numeric(as.character(ac..[ind2])))
        df <- df[rowSums(df[,3:4])>0,]
        df
    }

    getHet <- function(snp1,snp2)
    {
        gr1 <- GRanges(gsub("chr","",snp1[,1]),
                       IRanges(snp1[,2],snp1[,2]))
        gr2 <- GRanges(gsub("chr","",snp2[,1]),
                       IRanges(snp2[,2],snp2[,2]))
        ovs <- findOverlaps(gr1,gr2)
        snp1[queryHits(ovs),]
    }

    getAC <- function(ac_counts_paths, phases)
    {
        ## ##############################
        acs <- do.call("rbind",lapply(ac_counts_paths,function(dd)
        {
            gc()
            ok <- getHet(as.data.frame(data.table::fread(dd)),phases)
        }))
    }

    fitBinom.1dist <- function(counts, depths, steps=NULL, maxdepth=1000)
    {
        if(is.null(steps))
            steps <- if(length(counts)%/%3>10) 5 else 3
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        haploblocks <- cut(1:length(counts),pmax(length(counts)%/%steps,2))
        if(length(counts)<10) haploblocks <- rep(1,length(counts))
        lcounts <- split(counts,haploblocks)
        ldepths <- split(depths,haploblocks)
        lcounts <- lapply(1:length(lcounts), function(x) if(rnorm(1)<0) ldepths[[x]]-lcounts[[x]] else lcounts[[x]])
        counts <- sapply(lcounts,sum)
        depths <- sapply(ldepths,sum)
        counts[depths>maxdepth] <- round(counts[depths>maxdepth]/depths[depths>maxdepth]*maxdepth)
        depths[depths>maxdepth] <- maxdepth
        values <- seq(.5,1,0.001)
        llh <- sapply(values,function(x)
        {
            sum(log(dbinom(counts,size=depths,prob=x, log=F)+dbinom(counts,size=depths,prob=1-x, log=F)))
        })
        llh <- llh-max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        ##baf <- sum(values*normalised)
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs>.05)[1],which(cs>.95)[1])]
        c(q95[1],baf,q95[2])
    }

    fitBinom.1dist.noswitch <- function(counts, depths)
    {
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        values <- seq(0,1,0.001)
        llh <- sapply(values,function(x)
        {
            sum(dbinom(counts,size=depths,prob=x, log=T))
        })
        llh <- llh-max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        ##baf <- sum(values*normalised)
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs>.05)[1],which(cs>.95)[1])]
        c(q95[1],baf,q95[2])
    }
    is_distant_enough <- function(positions, distance=1000)
    {
        ##ord <- order(positions,decreasing=F)
        ##psort <- positions[ord]
        ##keep <- c(T,diff(psort)>distance)
        ##return(keep[order(ord,decreasing=F)])
        spositions <- sort(positions)
        keep <- numeric(0)
        prev_kept <- -Inf
        for (pos in spositions)
        {
            if (pos - prev_kept >= distance)
            {
                keep <- c(keep, pos)
                prev_kept <- pos
            }
        }
        positions%in%keep
    }
    getProfile <- function(df,
                           prof,
                           purity,
                           ploidy,
                           purs,
                           ploidies,
                           steps=NULL)
    {
        nprof <- data.frame(chr=as.character(prof[,"chromosome"]),
                            startpos=as.numeric(prof[,"start"]),
                            endpos=as.numeric(prof[,"end"]),
                            total_copy_number=as.numeric(prof[,"total_copy_number"]),
                            total_copy_number_logr=as.numeric(prof[,"total_copy_number_logr"]),
                            logr=as.numeric(prof[,"logr"]),
                            logr.sd=as.numeric(prof[,"logr.sd"]),
                            fitted=as.numeric(prof[,"total_copy_number"]),
                            q05=as.numeric(NA),
                            BAF=as.numeric(NA),
                            q95=as.numeric(NA),
                            q05_noswitch=as.numeric(NA),
                            BAF_noswitch=as.numeric(NA),
                            q95_noswitch=as.numeric(NA),
                            stringsAsFactors=F)
        nprof[nprof[,1]=="chr23",1] <- "chrX"
        nprof[nprof[,1]=="23",1] <- "X"
        nprof[nprof[,1]=="chr24",1] <- "chrY"
        nprof[nprof[,1]=="24",1] <- "Y"
        grseg <- GRanges(gsub("chr","",nprof[,"chr"]),
                         IRanges(as.numeric(nprof[,"startpos"]),
                                 as.numeric(nprof[,"endpos"])))
        grsnp <- GRanges(gsub("chr","",df[,1]),
                         IRanges(df[,2],
                                 df[,2]))
        ovs <- findOverlaps(grseg,
                            grsnp)
        qH <- queryHits(ovs)
        sH <- subjectHits(ovs)
        df[,3] <- as.numeric(as.character(df[,3]))
        df[,4] <- as.numeric(as.character(df[,4]))
        for(i in unique(qH))
        {
            inds <- 1:nrow(df)%in%sH[qH==i]
            inds[inds] <- inds[inds] & is_distant_enough(df[inds,2])
            if(sum(inds)>1)
            {
                nprof[i,c("q05","BAF","q95")] <- fitBinom.1dist(df[inds, 3],
                                                                rowSums(df[inds, c(3,4)]),
                                                                steps=steps)
                nprof[i,c("q05_noswitch","BAF_noswitch","q95_noswitch")] <- fitBinom.1dist.noswitch(df[inds, 3],
                                                                                                    rowSums(df[inds, c(3,4)]))
            }
        }
        nona <- !is.na(nprof[,"BAF"])
        sG <- .searchGrid(as.numeric(nprof[nona,"BAF"]),
                         as.numeric(nprof[nona,"logr"]),
                         as.numeric(nprof[nona,"endpos"])-as.numeric(nprof[nona,"startpos"]),
                         purs=purs,
                         ploidies=ploidies)
        sG.fixed <- .searchGrid(nprof[nona,"BAF"],
                               nprof[nona,"logr"],
                               nprof[nona,"endpos"]-nprof[nona,"startpos"],
                               purs=purity,
                               ploidies=ploidy)
        exceptNA <- function(vec,nona)
        {
            newvec <- rep(NA,length(nona))
            newvec[nona] <- vec
            newvec
        }
        nprof[,"ntot_free"] <- transform_bulk2tumour(nprof[,"logr"],sG$purity,sG$ploidy)
        nprof[,"ntot_fixed"] <- transform_bulk2tumour(nprof[,"logr"],sG.fixed$purity,sG.fixed$ploidy)
        list(nprof.free=cbind(nprof,
                              nA=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity),
                              nB=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity, retNa=FALSE),
                              nA_sc_raw=exceptNA(sG$profile$Nb,nona),
                              nB_sc_raw=exceptNA(sG$profile$Na,nona),
                              nA_sc=exceptNA(round(sG$profile$Nb),nona),
                              nB_sc=exceptNA(round(sG$profile$Na),nona)),
             purity.free=sG$purity,
             ploidy.free=sG$ploidy,
             nprof.fixed=cbind(nprof,
                               nA=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity),
                               nB=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity, retNa=FALSE),
                               nA_sc_raw=exceptNA(sG.fixed$profile$Nb,nona),
                               nB_sc_raw=exceptNA(sG.fixed$profile$Na,nona),
                               nA_sc=exceptNA(round(sG.fixed$profile$Nb),nona),
                               nB_sc=exceptNA(round(sG.fixed$profile$Na),nona)),
             purity.fixed=sG.fixed$purity,
             ploidy.fixed=sG.fixed$ploidy)
        }

    getAS_CNA_sample <- function(track,
                                 profile,
                                 ac_counts_paths,
                                 purs,
                                 ploidies,
                                 purity,
                                 ploidy,
                                 phases=NULL,
                                 path_to_phases=NULL,
                                 steps=NULL)
    {
        if(is.null(phases))
            phases <- readPhases(path_to_phases)
        ac <- getAC(ac_counts_paths, phases)
        ac.ph <- getPhasedInfo(ac, phases)
        prof <- getProfile(ac.ph,
                           prof=profile,
                           steps=steps,
                           purity=purity,
                           ploidy=ploidy,
                           purs=purs,
                           ploidies=ploidies)
        prof
    }


    phases <- NULL
    if(length(path_to_phases)==1)
    {
        print("## read-in Phases")
        phases <- readPhases(path_to_phases[[1]])
    }
    print("## derive Allele-specific Profiles")
    res$allProfiles_AS <- parallel::mclapply(1:length(res$allTracks.processed), function(x)
    {
        cat(".")
        getAS_CNA_sample(track=res$allTracks.processed[[x]],
                         profile=res$allProfiles[[x]],
                         ac_counts_paths=list_ac_counts_paths[[x]],
                         phases=phases,
                         purity=if(any(grepl("refitted",names(res)))) res$allSolutions.refitted.auto[[x]]$purity
                                else res$allSolutions[[x]]$purity,
                         ploidy=if(any(grepl("refitted",names(res)))) res$allSolutions.refitted.auto[[x]]$ploidy
                                else res$allSolutions[[x]]$ploidy,
                         purs=purs[[x]],
                         ploidies=ploidies[[x]],
                         path_to_phases=if(length(path_to_phases)>1) path_to_phases[[x]] else NULL,
                         steps=steps)
    },mc.cores=mc.cores)
    print("## write to disk and plot Allele-specific Profiles")
    pdf(paste0(outdir,"/all_as_cna_profiles_",projectname,".pdf"),width=15,height=5)
    tnull <- lapply(1:length(res$allProfiles_AS), function(x)
    {
        try({
            plot_AS_profile(res$allProfiles_AS[[x]]$nprof.fixed)
            title(paste0(names(res$allTracks)[x]," - bam",x) ,cex=.5)
        })
        write.table(res$allProfiles_AS[[x]]$nprof.fixed,
                    sep="\t",col.names=T,row.names=F,quote=F,
                    file=paste0(outdir,"/as_cna_profile_",names(res$allTracks)[x],"_bam",x,".txt"))
    })
    dev.off()
    res
}
