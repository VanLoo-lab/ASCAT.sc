sc_getBreakpoints_multipcf <- function(tmpdata, svinput, penalties, normalize)
{
    rmcc <- which(colnames(tmpdata)%in%c("start","end"))
    break_data <- function(tmpdata, sv, rmcc)
    {
        grcn <- GRanges(tmpdata[,"chr"],IRanges(tmpdata[,"start"],tmpdata[,"end"]))
        bpsv <- sort(sv[,2], decreasing=F)
        grregions <- GRanges(rep(tmpdata[1,"chr"],nrow(sv)+1),IRanges(c(0,bpsv),c(bpsv-1,1000000000)))
        ovs <- findOverlaps(grcn, grregions, select = "first")
        indices <- lapply(1:length(grregions),function(x)
        {
            which(ovs==x)
        })
        indices <- lapply(which(sapply(indices,length)>0),function(x) indices[[x]])
        ltmpdata <- lapply(indices, function(x) tmpdata[x,-rmcc])
    }
    ltmpdata <- break_data(tmpdata, svinput, rmcc)
    out <- lapply(penalties,function(penalty)
    {
        ok <- do.call("rbind",lapply(ltmpdata, function(tmpdata)
        {
            invisible(capture.output(suppressMessages({ok <- copynumber::multipcf(data=tmpdata,
                                                                                  gamma=penalty,
                                                                                  normalize=normalize)})))
            ok
        }))
        ok
    })
    out
}
