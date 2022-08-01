getlSe_excluded <- function(lSe, lInds)
{
    nlSe <- lapply(1:length(lSe),function(x)
    {
        if(length(lInds[[x]])>0)
            return(list(starts=lSe[[x]]$starts[-c(lInds[[x]])],
                        ends=lSe[[x]]$ends[-c(lInds[[x]])]))
        lSe[[x]]
    })
    names(nlSe) <- names(lSe)
    nlSe
}
