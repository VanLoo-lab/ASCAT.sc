getnlCTS_excluded <- function(nlCTS.tumour, lInds)
{
    nnlCTS.tumour <- lapply(1:length(nlCTS.tumour),function(chr)
    {
        if(length(lInds[[chr]])>0)
            return(nlCTS.tumour[[chr]][-lInds[[chr]],])
        nlCTS.tumour[[chr]]
    })
    names(nnlCTS.tumour) <- names(nlCTS.tumour)
    nnlCTS.tumour
}

