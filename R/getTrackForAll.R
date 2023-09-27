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
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS=lCTS,lSegs=lSegs)
    return(tracks)
}
