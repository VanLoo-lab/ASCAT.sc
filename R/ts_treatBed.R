ts_treatBed <- function(bed, add)
{
    bed <- bed[,1:3]
    newbed <- NULL
    chrs <- unique(bed[,1])
    grbed <- GRanges(bed[,1],IRanges(bed[,2]-add,bed[,3]+add))
    covs <- coverage(grbed)
    nms <- names(covs)
    covs <- covs>0
    cov1  <-  IRanges::slice(covs, lower=1)
    ngrbed  <-  as(cov1, "GRanges")
    newbed <- cbind(as.vector(seqnames(ngrbed)),start(ranges(ngrbed)),end(ranges(ngrbed)))
    newbed <- data.frame("chr"=newbed[,1],
                         "start"=as.numeric(newbed[,2]),
                         "end"=as.numeric(newbed[,3]))
    colnames(newbed) <- colnames(bed)
    newbed
}

