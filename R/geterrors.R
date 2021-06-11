geterrors <- function(rho,
                      phi,
                      meansSeg,
                      weights,
                      sds=1,
                      gamma=1,
                      ismale=F,
                      isPON=F,
                      isX=F,
                      distance="mse")
{
    if(distance[1]=="mse")
    {
        signal <- transform_bulk2tumour(meansSeg,rho,phi, gamma=gamma, isPON=isPON, ismale=ismale, isX=isX)
        return(mean(((round(signal)-signal)/sds)^2*weights/1000,na.rm=T))
    }
    if(grepl("stat",distance[1]))
    {
        ## experimental
        ## signal <- transform_bulk2tumour(meansSeg,rho,phi, gamma=gamma)
        ## return(mean(((round(signal)-signal)/sds)^2*weights/1000,na.rm=T))
        ## return(sum(dnorm(signal,mean=round(signal),sd=sds,log=T)*weights/1000,na.rm=T))
    }
}
