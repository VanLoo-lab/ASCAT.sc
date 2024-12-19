getTrackForAll <- function(bamfile,
                           window,
                           lCT=NULL,
                           lSe=NULL,
                           lGCT=NULL,
                           lNormals=NULL,
                           allchr=1:22,
                           sdNormalise=0,
                           segmentation_alpha=0.01,
			   isDuplicate=F,
			   isSecondaryAlignment=F,
			   isNotPassingQualityControls=NA,
			   isUnmappedQuery=NA,
			   mapqFilter=30,
                           doSmooth=TRUE,
                           SBDRY=NULL,
                           svinput=NULL,
                           doSeg=TRUE)
{
    if(is.null(lSe))
    {
        print("   ## get Start-End of segments")
        lSe <- lapply(allchr,function(chr) getStartsEnds(window,paste0(chr)))
    }
    ## ##################################################
    if(is.null(lCT))
    {
        print("   ## get Coverage Track")
        lCT <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                             chr=paste0(chr),
                                                             lSe[[chr]]$starts,
                                                             lSe[[chr]]$ends,
                                                             isDuplicate=isDuplicate,
                                                             isSecondaryAlignment=isSecondaryAlignment,
                                                             isNotPassingQualityControls=isNotPassingQualityControls,
                                                             isUnmappedQuery=isUnmappedQuery,
                                                             mapqFilter=mapqFilter))
    }
    if(is.null(lGCT))
    {
        print("   ## get GC content")
        lGCT <- lapply(allchr,function(chr) gcTrack(chr,lSe[[chr]]$starts,lSe[[chr]]$ends))
    }
    ## ##################################################
    if(doSmooth)
    {
        print("   ## correct for GC content")
        lCTS <- smoothCoverageTrack(lCT=lCT,lSe=lSe,lGCT=lGCT, lNormals=lNormals)
        names(lCTS) <- allchr
        gc()
    }
    ## ##################################################
    print("   ## segment Tracks")
    require(DNAcopy)
    if(is.null(svinput))
    {
        lSegs <- lapply(1:length(lCTS),function(x)
        {
            segments<- segmentTrack(lCTS[[x]]$smoothed,
                                    chr=paste0(x),
                                    sd=sdNormalise,
                                    lSe[[x]]$starts,
                                    lSe[[x]]$ends,
                                    SBDRY=SBDRY,
                                    ALPHA=segmentation_alpha)
        })
    }
    if(!is.null(svinput))
    {
        segmentTrack_by_parts<-function(smoothed,
                                        chr,
                                        allchr,
                                        sd,
                                        starts,
                                        ends,
                                        SBDRY,
                                        svinput,
                                        ALPHA)
        {
            widths <- ends-starts
            grcn <- GRanges(allchr[chr],IRanges(starts+round(widths/4),ends-round(widths/4)))
            bpsv <- sort(svinput[,2], decreasing=F)
            grregions <- GRanges(rep(allchr[chr],nrow(svinput)+1),IRanges(c(0,bpsv),c(bpsv-1,1000000000)))
            ovs <- findOverlaps(grcn, grregions, select = "first")
            indices <- lapply(1:length(grregions),function(x)
            {
                which(ovs==x)
            })
            indices <- lapply(which(sapply(indices,length)>0),function(x) indices[[x]])
            ll <- lapply(indices,function(ii)
            {
                segmentTrack(smoothed[ii],
                             chr=chr,
                             sd=sd,
                             starts[ii],
                             ends[ii],
                             SBDRY=SBDRY,
                             ALPHA=ALPHA)
            })
            out <- list(data = do.call("rbind",lapply(ll,function(x) x$data)),
                        output = do.call("rbind",lapply(ll,function(x) x$output)))
        }
        lSegs <- lapply(1:length(lCTS),function(x)
        {
            sv_condchr <- gsub("chr","",svinput[,1])==gsub("chr","",allchr[x])
            svinput. <- svinput[sv_condchr,]
            if(sum(sv_condchr)>0)
            {
                segments <- segmentTrack_by_parts(smoothed=lCTS[[x]]$smoothed,
                                                  chr=x,
                                                  allchr=allchr,
                                                  sd=sdNormalise,
                                                  starts=lSe[[x]]$starts,
                                                  ends=lSe[[x]]$ends,
                                                  SBDRY=SBDRY,
                                                  svinput=svinput.,
                                                  ALPHA=segmentation_alpha)
            }
            else
                segments <- segmentTrack(lCTS[[x]]$smoothed,
                                         chr=paste0(x),
                                         sd=sdNormalise,
                                         lSe[[x]]$starts,
                                         lSe[[x]]$ends,
                                         SBDRY=SBDRY,
                                         ALPHA=segmentation_alpha)
        })
    }
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS=lCTS,lSegs=lSegs)
    return(tracks)
}
