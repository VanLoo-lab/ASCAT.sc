sc_getMatrixOfCalls <- function(listTracks, lSolutions, lSe)
{
    allcalls <- lapply(1:length(listTracks),function(x) sc_getCalls(listTracks[[x]],
                                                                    lSolutions[[x]]$purity,
                                                                    lSolutions[[x]]$ploidy))
    mat <- lapply(allcalls,function(x) unlist(lapply(x,function(y)
    {
        v <- rep(y$CNtot,y$num.mark)
        v
    })))
    isNull <- sapply(mat,is.null)
    mat <- sapply(which(!isNull),function(x) mat[[x]])
    rnms <- unlist(lapply(1:length(allcalls[[1]]),function(x) paste0(x,":",lSe[[x]]$starts,"-",lSe[[x]]$ends)))
    rownames(mat) <- rnms
    mat
}
