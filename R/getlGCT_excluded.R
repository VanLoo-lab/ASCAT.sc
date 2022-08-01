getlGCT_excluded <- function(lGCT, lInds)
{
    nlGCT <- lapply(1:length(lGCT),function(x)
    {
        if(length(lInds[[x]])>0)
            return(lGCT[[x]][-c(lInds[[x]])])
        lGCT[[x]]
    })
    names(nlGCT) <- names(lGCT)
    nlGCT
}
