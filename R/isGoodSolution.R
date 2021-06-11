isGoodSolution <-
function(meansSeg,maxSize0=20)
{
    starts <- unlist(lapply(meansSeg,function(x) lapply(x,function(y) y$start)))
    ends <- unlist(lapply(meansSeg,function(x) lapply(x,function(y) y$end)))
    sizes <- (ends-starts)/1000000
    cn <- unlist(lapply(meansSeg,function(x) lapply(x,function(y) y$roundmu)))
    cn0 <- round(cn)<=0
    totalsizes0 <- sum(sizes[cn0],na.rm=T)
    if(totalsizes0<maxSize0) return(TRUE)
    FALSE
}
