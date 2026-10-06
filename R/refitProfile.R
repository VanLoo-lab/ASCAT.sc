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
                         CHRS=NULL,
                         gridpur=seq(-.05,.05,.01),
                         gridpl=seq(-.1,.2,.01))
{
    if(is.null(CHRS)) CHRS <- 1:22
    errs_orig <- solution$errs
    profile <- getProfile(fitProfile(track,solution$purity,solution$ploidy, gamma=gamma, ismale=ismale, isPON=isPON),
                          CHRS=CHRS)
    if(!is.na(chr1))
    {
        profile1 <- profile[profile[,"chromosome"]==chr1,,drop=F]
        if (nrow(profile1) == 0)
            stop(paste0("Chromosome '", chr1, "' has no segments in this profile -- choose a different chromosome"))
        # Use the longest segment that has data (non-NA logr): with sparse
        # single-cell coverage the longest segment by span can have no
        # covered bins.
        valid1 <- !is.na(profile1[,"logr"])
        if (!any(valid1))
            stop(paste0("Chromosome '", chr1, "' has no segments with usable data (all NA logR) -- try a different chromosome"))
        cand1 <- which(valid1)
        longestsegment1 <- cand1[which.max(profile1[cand1,"end"]-profile1[cand1,"start"])[1]]
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
        if (nrow(profile2) == 0)
            stop(paste0("Chromosome '", chr2, "' has no segments in this profile -- choose a different chromosome"))
        valid2 <- !is.na(profile2[,"logr"])
        if (!any(valid2))
            stop(paste0("Chromosome '", chr2, "' has no segments with usable data (all NA logR) -- try a different chromosome"))
        cand2 <- which(valid2)
        longestsegment2 <- cand2[which.max(profile2[cand2,"end"]-profile2[cand2,"start"])[1]]
    }
    if(is.na(chr2))
    {
        if(is.na(ind2)) stop("Specify either chromosome or index for second segment")
        profile2 <- profile
        longestsegment2 <- ind2
    }
    logr1 <- 2^(profile1[longestsegment1,"logr"]/gamma)
    logr2 <- 2^(profile2[longestsegment2,"logr"]/gamma)
    # A segment selected by index (ind1/ind2) can still have no data.
    if (is.na(logr1) || is.na(logr2))
        stop("Selected segment has no underlying data (NA logR) -- try a different chromosome, or one with more covered bins")
    purity <- (2*logr1/total1/logr2-2/total1)/(1-logr1*total2/logr2/total1+logr1/logr2/total1*2-2/total1)
    ploidy <- (total2*purity+(1-purity)*2)/logr2
    gridpur <- purity+gridpur
    gridpl <- ploidy+gridpl
    # is.finite() also drops NA/NaN/Inf, which a range check alone
    # would keep.
    gridpur <- gridpur[is.finite(gridpur) & gridpur>0 & gridpur<=1]
    gridpl  <- gridpl[is.finite(gridpl) & gridpl>0]
    if(length(gridpur)==0 | length(gridpl)==0) stop("Not possible: ploidy<0 or purity \u2209 [0,1]")
    newsol <- searchGrid(track,
               purs=gridpur,
               ploidies=gridpl,
               gamma=gamma,
               ismale=ismale,
               isPON=isPON)
    newsol$errs <- errs_orig
    newsol
}
