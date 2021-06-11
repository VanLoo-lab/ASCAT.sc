removeCentromeres <- function(lSe, centromeres, chr)
{
    require(GenomicRanges)
    if(is.null(centromeres)) return(lSe)
    grlse <- GRanges(rep(chr,length(lSe$starts)),
                     IRanges(lSe$starts,
                             lSe$ends))
    ov <- findOverlaps(grlse,centromeres)
    list(starts=lSe$starts[-c(queryHits(ov))],
         ends=lSe$ends[-c(queryHits(ov))])
}
