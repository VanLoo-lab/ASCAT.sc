treatlSe <- function(lSe, window=ceiling(as.numeric(WINDOW)/as.numeric("10000")))
{
    nlSe <- lapply(lSe,function(x)
    {    
        print("beforenchar2")
        l <- nchar(x$starts)
        starts <- getstartends(end=l,window=window)$starts
        ends <- getstartends(end=l,window=window)$ends
        print("afternchar2")
        list(starts=x$starts[starts],
             ends=x$ends[ends])
    })
    names(nlSe) <- names(lSe)
    nlSe
}
