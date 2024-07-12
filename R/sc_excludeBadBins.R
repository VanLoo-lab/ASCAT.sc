sc_excludeBadBins <- function(res)
{
    if(!is.null(res$lCTS.normal$nlCTS.normal))
    {
        lInds <- filterBins(allTracks=res$lCTS.normal$nlCTS.normal, logr=NULL, lSe=res$nlSe, IQRC=1.5)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        res$nlCTS.tumour <- getnlCTS_excluded(res$nlCTS.tumour, lInds)
    }
    res
}
