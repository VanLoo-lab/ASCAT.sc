<<<<<<< HEAD
meth_winsorise_ascat <- function(x, musmooth=F, window=9)
=======
meth_winsorise_ascat <- function(x)
>>>>>>> d619bd95c4c64791b532ea98ba54d68877ef031b
{
    .madWins <- function(x,tau,k)
    {
        .medianFilter <- function(x,k){
            n <- length(x)
            filtWidth <- 2*k + 1
                                        #Make sure filtWidth does not exceed n
            if(filtWidth > n){
                if(n==0){
                    filtWidth <- 1
                }else if(n%%2 == 0){
                                        #runmed requires filtWidth to be odd, ensure this:
                    filtWidth <- n - 1
                }else{
                    filtWidth <- n
                }
            }
            runMedian <- stats::runmed(x,k=filtWidth,endrule="median")
            return(runMedian)
        }
        .psi <- function(xwin, zplus, zminus)
        {
            xwin[xwin < -zminus] <- NA
            xwin[xwin > zplus] <- NA
            return(xwin)
        }
        .winsorise <- function(d, xhat)
        {
            SDplus <- median(d[d>=0])
            SDminus <- median(abs(d[d<=0]))
            zplus <- tau*SDplus
            zminus <- tau*SDminus
            xwin <- xhat + .psi(d, zplus, zminus)
        }
        xhat <- .medianFilter(x,k)
        d <- x-xhat
        qq <- quantile(xhat,prob=seq(0,1,.05))[-c(1)]
        bools <- rep(F,length(xhat))
        xwin <- d
        for(q in qq)
        {
            bools <- xor(bools,xhat<=q)
            xwin[bools] <- .winsorise(d[bools], xhat[bools])
        }
        return(list(ywin=xwin))
    }
<<<<<<< HEAD
    if(musmooth) x <- caTools::runmean(x,k=window)
    lrwins = vector(mode="numeric",length=length(x))
    lrwins[is.na(x)] = NA
    lrwins[!is.na(x)] = .madWins(x[!is.na(x)],4,50)$ywin
=======
    lrwins = vector(mode="numeric",length=length(x))
    lrwins[is.na(x)] = NA
    lrwins[!is.na(x)] = .madWins(x[!is.na(x)],6,50)$ywin
>>>>>>> d619bd95c4c64791b532ea98ba54d68877ef031b
    lrwins
}
