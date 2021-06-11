segmentTrack <- function (covtrack,
                          chr,
                          starts,
                          ends = NA,
                          sd = 0,
                          min.width = 5,
                          ALPHA=0.01,
                          smooth=T)
{
    require(DNAcopy)
    covtrack <- covtrack * rnorm(length(covtrack), mean = 1,
                                 sd = sd)
    cna <- CNA(covtrack, chr = rep(chr, length(covtrack)), maploc = starts,
               data.type = "logratio")
    if(smooth)
        cna <- smooth.CNA(cna)
    capture.output(cna <- segment(cna, min.width = min.width,alpha=ALPHA))
    cna
}
