transform_bulk2tumour <- function(x, rho, psi, ismale=F, isX=F, isPON=F, gamma=1)
{
    isdivideby2 <- as.numeric(ismale)*as.numeric(isX)*as.numeric(isPON)
    divideby2 <- 1/(isdivideby2+1)
    baseline <- ifelse(ismale & isX, 1, 2)
    return((psi*2^(x/gamma)-baseline*(1-rho))/rho*divideby2)
}
