meth_smoothAllChr <- function(logr,starts,ends,chrs,step=5)
{
    nl <- ns <- ne <- nc <- NULL
    for(i in unique(chrs))
    {
        keepchr <- chrs==i
        ss <- meth_smooth(logr[keep],starts[keepchr],ends[keepchr],step=step)
        nl <- c(nl,ss$nlg)
        ns <- c(ns,ss$starts)
        ne <- c(ne,ss$ends)
        nc <- c(nc,rep(i,length(ss$nlg)))
    }
    list(nl,ns,ne,nc)
}
