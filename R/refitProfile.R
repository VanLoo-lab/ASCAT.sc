refitProfile <- function(track,
                         solution,
                         chr1=NA,
                         ind1=NA,
                         total1,
                         chr2=NA,
                         ind2=NA,
                         total2,
                         gamma=1,
                         ismale=F,
                         isPON=F,
                         gridpur=seq(-.05,.05,.01),
                         gridpl=seq(-.1,.2,.01))
{
    profile <- getProfile(fitProfile(track,solution$purity,solution$ploidy, gamma=gamma, ismale=ismale, isPON=isPON))
    if(!is.na(chr1))
    {
        profile1 <- profile[profile[,"chromosome"]==chr1,,drop=F]
        longestsegment1 <- which.max(profile1[,"end"]-profile1[,"start"])[1]
    }
    if(is.na(chr1))
    {
        if(is.na(ind1)) stop("Specify either chromosome or index for first segment")
        profile1 <- profile
        longestsegment1 <- ind1
    }
    if(!is.na(chr2))
    {
        profile2 <- profile[profile[,"chromosome"]==chr2,,drop=F]
        longestsegment2 <- which.max(profile2[,"end"]-profile2[,"start"])[1]
    }
    if(is.na(chr2))
    {
        if(is.na(ind2)) stop("Specify either chromosome or index for second segment")
        profile2 <- profile
        longestsegment2 <- ind2
    }
    logr1 <- 2^(profile1[longestsegment1,"logr"]/gamma)
    logr2 <- 2^(profile2[longestsegment2,"logr"]/gamma)
    purity <- (2*logr1/total1/logr2-2/total1)/(1-logr1*total2/logr2/total1+logr1/logr2/total1*2-2/total1)
    ploidy <- (total2*purity+(1-purity)*2)/logr2
    gridpur <- purity+gridpur
    gridpl <- ploidy+gridpl
    gridpur <- gridpur[gridpur>0 & gridpur<=1]
    gridpl <- gridpl[gridpl>0]
    if(length(gridpur)==0 | length(gridpl)==0) stop("Not possible: ploidy<0 or purity ∉ [0,1]")
    searchGrid(track,
               purs=gridpur,
               ploidies=gridpl,
               gamma=gamma,
               ismale=ismale,
               isPON=isPON)
}
