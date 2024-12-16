sc_excludeBadBins <- function(res)
{
    if(!is.null(res$lCTS.normal$nlCTS.normal))
    {
        print("Using normal samples for removal of bad bins")
        lInds <- filterBins(allTracks=res$lCTS.normal$nlCTS.normal, logr=NULL, lSe=res$nlSe, IQRC=1.5)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        if(!is.null(res$lNormals))
            res$lNormals <- lapply(names(res$lNormals), function(x) getnlCTS_excluded(res$lNormals[[x]], lInds))
        res$allTracks <- lapply(names(res$allTracks),function(x)
        {
            res$allTracks[[x]]$nlCTS.tumour = getnlCTS_excluded(res$allTracks[[x]]$nlCTS.tumour, lInds)
            res$allTracks[[x]]
        })
    }
    else
    {
        print("Using all samples for removal of bad bins")
        allTracks <- lapply(res$allTracks, function(x) {
            list(lCTS = x$nlCTS.tumour)})
        lInds <- filterBins(allTracks=allTracks, logr=NULL, lSe=res$nlSe, fraction=.1, k=11)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        if(!is.null(res$lNormals))
            res$lNormals <- lapply(names(res$lNormals), function(x) getnlCTS_excluded(res$lNormals[[x]], lInds))
        res$allTracks <- lapply(names(res$allTracks),function(x)
        {
            res$allTracks[[x]]$nlCTS.tumour = getnlCTS_excluded(res$allTracks[[x]]$nlCTS.tumour, lInds)
            res$allTracks[[x]]
        })
    }
    res
}
