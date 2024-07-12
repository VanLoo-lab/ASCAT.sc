meth_winsorise_ascat <- function(x, apply_wins=FALSE)
{
    if(!apply_wins) return(x)
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
            xwin[xwin < -zminus] <- zminus
            xwin[xwin > zplus] <- zplus
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
    lrwins = vector(mode="numeric",length=length(x))
    lrwins[is.na(x)] = NA
    lrwins[!is.na(x)] = .madWins(x[!is.na(x)],4,50)$ywin
    lrwins
}
