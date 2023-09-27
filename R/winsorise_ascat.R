winsorise_ascat <- function(x)
{
    madWins <- function(x,tau,k)
    {
        medianFilter <- function(x,k){
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
        xhat <- medianFilter(x,k)
        d <- x-xhat
        SD <- mad(d)
        z <- tau*SD
        xwin <- xhat + psi(d, z)
        outliers <- rep(0, length(x))
        outliers[x > xwin] <- 1
        outliers[x < xwin] <- -1
        return(list(ywin=xwin,sdev=SD,outliers=outliers))
    }
    psi <- function(x,z)
    {
        xwin <- x
        xwin[x < -z] <- -z
        xwin[x > z] <- z
        return(xwin)
    }
    lrwins = vector(mode="numeric",length=length(x))
    lrwins[is.na(x)] = NA
    lrwins[!is.na(x)] = madWins(x[!is.na(x)],3,15)$ywin
    lrwins
}

