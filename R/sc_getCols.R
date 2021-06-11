sc_getCols <-
function(mmm,channel="R")
{
    poss <- as.character(c(-2,-1,0,1,2))
    Rs <- c(.6,.3,0,0,0)
    Gs <- c(0,0,0,.3,.6)
    Bs <- c(.2,.2,0,.5,.5)
    names(Rs) <- names(Gs) <- names(Bs) <- poss
    get(paste0(channel,"s"))[as.character(mmm)]
}
