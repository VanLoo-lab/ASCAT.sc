getAS_CNA_smoothed <- function(res,
                               mc.cores)
{
    fitIntegers.2D <- function(baf, ntot, iter=1, INDEX)
    {
        NNN=length(baf)
        baf=c(baf,1-baf)
        ntot=c(ntot,ntot)
        assignRemainingT <- function(ntot, remaining, combi, gaussian.mean, gaussian.sd)
        {
            print(paste0("warning: should not be here - ", min(sum(remaining),length(remaining))))
            ntot <- ntot[remaining]
            pgroup <- matrix(0,length(ntot),nrow(combi))
            for(i in 1:nrow(combi))
            {
                pgroup[,i] <- dnorm(ntot, gaussian.mean[[i]], gaussian.sd[[i]],log=T)
            }
            apply(pgroup,1,which.max)
        }
        assignRemaining <- function(ntot, baf, remaining, combi, gaussian.mean, gaussian.sd)
        {
            assT <- assignRemainingT(ntot, remaining, combi, gaussian.mean, gaussian.sd)
            ntot <- ntot[remaining]
            baf <- baf[remaining]
            pgroup <- matrix(0,length(ntot),nrow(combi))
            for(i in 1:nrow(combi))
            {
                pgroup[,i] <- (ntot-gaussian.mean[[i]])^2 + (baf-combi[i,1]/(combi[i,1]+combi[i,2]))^2
            }
            mywhich.min <- function(...)
            {
                ret <- which.min(...)
                if(length(ret)==0) return(NA)
                ret
            }
            ww <- apply(pgroup,1,mywhich.min)
            print(sum(is.na(ww)))
            if(sum(is.na(ww))>0) ww[is.na(ww)] <- assT[is.na(ww)]
            ww
        }
        ## ############
        AS_mode_on <- TRUE
        if(mean(!is.na(baf))<.5)
        {
            baf <- rep(0.001,length(baf))
            AS_mode_on <- FALSE
        }
        rmNA <- !is.na(baf) & !is.na(ntot)
        BAF <- baf
        NTOT <- ntot
        ## ############
        baf <- baf[rmNA]
        ntot <- ntot[rmNA]
        baf[baf==0] <- 0.001
        baf[baf==1] <- 1-0.001
        allele1 <- pmax(0,round(baf*round(ntot)))
        allele2 <- pmax(0,round(ntot)-allele1)
        ## ############
        combi <- unique(cbind(allele1,allele2))
        beta.sd <- lapply(1:nrow(combi),function(x) NULL)
        gaussian.sd <- lapply(1:nrow(combi),function(x) NULL)
        for(j in 1:iter)
        {
            ##cat(".")
            beta.shape1 <- list()
            beta.shape2 <- list()
            beta.mean1 <- list()
            gaussian.mean <- list()
            pgroup <- matrix(0,length(baf),nrow(combi))
            for(i in 1:nrow(combi))
            {
                ww <- allele1==combi[i,1] & allele2==combi[i,2]
                beta.mean <- combi[i,1]/(combi[i,1]+combi[i,2])
                beta.mean[beta.mean>.99] <- .99
                beta.mean[beta.mean<.01] <- .01
                beta.sd. <- min(max(.0015,var(baf[ww]),na.rm=T),.015*(1-abs(beta.mean-0.5)),na.rm=T)
                beta.sd[[i]] <- if(is.null(beta.sd[[i]])) beta.sd. else if(beta.sd[[i]]>beta.sd.) beta.sd[[i]] else beta.sd.
                beta.shape1[[i]] <- ((1-beta.mean)/beta.sd[[i]]-1/beta.mean)*beta.mean^2
                beta.shape2[[i]] <- beta.shape1[[i]]*(1/beta.mean-1)
                gaussian.mean[[i]] <- sum(combi[i,c(1,2)])
                ##gaussian.sd[[i]] <- min(max(0.05,sd(ntot[ww])),0.25)
                ##gaussian.sd. <- min(max(0.05,sd(ntot[ww]),na.rm=T),0.125,na.rm=T)
                gaussian.sd. <- 0.5/3
                gaussian.sd[[i]] <- if(is.null(gaussian.sd[[i]])) gaussian.sd. else if(gaussian.sd[[i]]>gaussian.sd.) gaussian.sd[[i]] else gaussian.sd.
                ##gaussian.sd[[i]] <- sd(ntot[ww])
                pBAF <- normalize_llh(dbeta(baf, beta.shape1[[i]], beta.shape2[[i]],log=T))
                pLOGR <- normalize_llh(dnorm(ntot, gaussian.mean[[i]], gaussian.sd[[i]],log=T))
                pgroup[,i] <- pBAF+pLOGR*1.5
            }
            ww <- apply(pgroup,1,which.max)
            counts <- table(ww)[as.character(1:nrow(combi))]
            counts[is.na(counts)] <- 0
            allele1 <- combi[ww,1]
            allele2 <- combi[ww,2]
            if(j<iter)
            {
                combi <- combi[counts>0,,drop=F]
                if(nrow(combi)==1)
                {
                    combi <- rbind(allele1=round(median(baf)*round(median(ntot))),
                                   allele2=round(ntot)-round(median(baf)*round(median(ntot))))
                    allele1 <- round(BAF*round(NTOT))
                    allele2 <- round(NTOT)-allele1
                    return(list(ww=rep(1,length(NTOT)),
                                baf=BAF,
                                ntot=NTOT,
                                allele1Inferred=combi[1,1],
                                allele2Inferred=combi[1,2],
                                allele1=allele1,
                                allele2=allele2,AS_mode_on=AS_mode_on))
                }
                beta.sd <- lapply(which(counts>0),function(x) beta.sd[[x]])
                gaussian.sd <- lapply(which(counts>0),function(x) gaussian.sd[[x]])
            }
        }
        ## ############
        ww[apply(pgroup,1,function(x)
        {
            ord <- sort(x,decreasing=T)[1:2]
            all(is.infinite(x) & ord[1]==ord[2])
        })] <- NA
        ## ############
        WW <- rep(NA,length(NTOT))
        WW[which(rmNA)] <- ww
        remaining <- which(is.na(WW))
        if(length(remaining)>0)
        {
            assRemain <- assignRemaining(NTOT, BAF, remaining, combi, gaussian.mean, gaussian.sd)
            if(length(assRemain)==0)
            {
                print(INDEX)
            }
            WW[remaining] <- assRemain
        }
        ## ############
        allele1 <- round(BAF*round(NTOT))
        allele2 <- round(NTOT)-allele1
        A1 <- tapply(1:length(BAF),WW,function(x) round(median(BAF[x],na.rm=T)*round(median(NTOT[x],na.rm=T))))
        allele1Inferred <- A1[as.character(WW)]
        A2 <- tapply(1:length(BAF),WW,function(x) round(median(round(NTOT[x])-allele1Inferred[x],na.rm=T)))
        allele2Inferred <- A2[as.character(WW)]
        list(ww=WW[1:NNN],
             baf=BAF[1:NNN],
             ntot=NTOT[1:NNN],
             allele1Inferred=allele1Inferred[1:NNN],
             allele2Inferred=allele2Inferred[1:NNN],
             allele1=allele1[1:NNN],
             allele2=allele2[1:NNN],
             AS_mode_on=AS_mode_on[1:NNN])
    }

    lProfs <- res$allProfiles_AS
    if(grepl("filters",names(res)))
    {
        nms <- names(lProfs)
        lProfs <- lapply(which(res$filters), function(x)
        {
            lProfs[[x]]
        })
        names(lProfs) <- nms[which(res$filters)]
    }
    alleles <- parallel::mclapply(1:nrow(lProfs[[1]]),function(index)
    {
        cat(".")
        baf <- sapply(lProfs,function(x) x[index,"BAF"])
        ntot <- sapply(lProfs,function(x) x[index,"ntot"])
        fitted <- fitIntegers.2D(baf, ntot, iter=100, INDEX=index)
    },mc.cores=mc.cores)
    newlProfs <- lapply(1:length(lProfs),function(x)
    {
        tt <- cbind(lProfs[[x]],
                    allele1Inferred=sapply(alleles,function(y) y$allele1Inferred[x]),
                    allele2Inferred=sapply(alleles,function(y) y$allele2Inferred[x]),
                    AS_mode_ON=sapply(alleles,function(y) y$AS_mode_on[1]))
    })
    names(newlProfs) <- names(lProfs)
    res$allProfiles_AS_smoothed <- newlProfs
    res
}
