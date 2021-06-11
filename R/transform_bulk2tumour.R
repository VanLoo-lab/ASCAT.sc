transform_bulk2tumour <- function(x,rho,phi, ismale=F, isX=F, isPON=F, gamma=1)
{
    isdivideby2 <- as.numeric(ismale)*as.numeric(isX)*as.numeric(isPON)
    divideby2 <- 1/(isdivideby2+1)
    return((phi*2^(x/gamma)-2*(1-rho))/rho*divideby2)
}
