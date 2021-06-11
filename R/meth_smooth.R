meth_smooth <- function(logr, starts, ends, step=20, is_isolength=T)
{
    isolength <- function(logr, starts, ends, step=20)
    {
        LL <- length(starts)
        nbins <- LL%/%step
        start.bins <- round(seq(starts[1],starts[LL],length.out=nbins))
        end.bins <- c(start.bins[2:nbins]-1,ends[LL])
        gr.bins <- GenomicRanges::GRanges(rep("chr1",nbins),IRanges(start.bins,end.bins))
        gr.probes <- GenomicRanges::GRanges(rep("chr1",LL),IRanges(starts,width=1))
        overlaps <- findOverlaps(gr.bins,gr.probes)
        keep <- 1:nbins%in%queryHits(overlaps)
        value.bins <- tapply(subjectHits(overlaps),queryHits(overlaps),function(x) median(logr[x],na.rm=T))
        value.bins <- value.bins[order(as.numeric(names(value.bins)),decreasing=F)]
        nprobes <- table(queryHits(overlaps))
        return(data.frame(nlg=as.numeric(value.bins),
                          starts=start.bins[keep],
                          ends=end.bins[keep],
                          nprobes=as.vector(nprobes[names(value.bins)])))
    }
    mergebins <- function(isl,minprobes=15,maxsize=5000000,startclip=1000000,endclip=1000000)
    {
        while(isl$nprobes[1]<minprobes)
        {
            isl$ends[1] <- isl$ends[2]
            isl$nlg[1] <- sum(isl$nlg[1:2]*isl$nprobes[1:2]/sum(isl$nprobes[1:2],na.rm=T),na.rm=T)
            isl$nprobes[1] <- sum(isl$nprobes[1:2])
            isl <- isl[-2,]
        }
        while(any(isl$nprobes<minprobes))
        {
            ww <- which.min(isl$nprobes)[1]
            isl$ends[ww-1] <- isl$ends[ww]
            isl$nlg[ww-1] <- sum(isl$nlg[(ww-1):ww]*isl$nprobes[(ww-1):ww]/sum(isl$nprobes[(ww-1):ww]),na.rm=T)
            isl$nprobes[ww-1] <- sum(isl$nprobes[(ww-1):ww])
            isl <- isl[-ww,]
        }
        isl <- isl[isl$ends-isl$starts<maxsize,]
        isl <- isl[isl$starts>min(isl$starts)+startclip & isl$ends<max(isl$ends)-endclip,]
        isl
    }
    isoprobes <- function(logr, starts, ends, step=20)
    {
        starti <- seq(1,length(starts)-step,step)
        endi <- c(starti[-c(1)]-1,length(starts))
        nlg <- sapply(1:length(starti),function(x) mean(logr[starti[x]:endi[x]],na.rm=T))
        return(list(nlg=nlg,
                    starts=starts[starti],
                    ends=ends[endi]))
    }
    if(is_isolength)
    {
        isl <- isolength(logr, starts, ends, step=20)
        misl <- mergebins(isl)
        return(misl)
    }
    else
    {
        return(isoprobes(logr, starts, ends, step=20))
    }
}
