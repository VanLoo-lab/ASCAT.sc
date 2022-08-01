getlInds <- function(lSe, lExclude)
{
    lInds <- lapply(1:length(lSe), function(x)
    {
        gr1 <- GRanges(x,IRanges(lSe[[x]]$starts,lSe[[x]]$ends))
        gr2 <- GRanges(x,IRanges(lExclude[[x]]$starts, lExclude[[x]]$ends))
        ovs <- findOverlaps(gr1,gr2)
        queryHits(ovs)
    })
}
