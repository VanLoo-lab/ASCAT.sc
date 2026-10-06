refitProfile_shift <- function (track,
                                solution,
                                gamma = 1,
                                ismale = F,
                                isPON = F,
                                CHRS = NULL,
                                shift = c(-1,1),
                                ismedian=FALSE,
                                gridpur = seq(-0.05,
                                              0.05, 0.01),
                                gridpl = seq(-0.1, 0.2, 0.01))
{
    shift <- shift[1]
    if(is.null(CHRS)) CHRS <- 1:22
    errs_orig <- solution$errs
    profile <- getProfile(fitProfile(track, solution$purity,
                                     ismedian=ismedian,
                                     solution$ploidy, gamma = gamma,
                                     ismale = ismale, isPON = isPON),
                          CHRS=CHRS)
    sizes <- (as.numeric(profile[,"end"])-as.numeric(profile[,"start"]))/1000000
    meanlogr <- tapply(1:nrow(profile),profile[,"total_copy_number"],function(x)
    {
        sum(sizes[x]*as.numeric(profile[x,"logr"]))/sum(sizes[x])
    })
    lengthlogr <- tapply(1:nrow(profile),profile[,"total_copy_number"],function(x) sum(sizes[x]))
    lengthlogr <- lengthlogr[names(meanlogr)]
    # Keep only copy-number states that stay above 0 once shifted.
    keep <- (as.numeric(names(lengthlogr)) + shift) > 0
    lengthlogr <- lengthlogr[keep]
    meanlogr <- meanlogr[keep]
    # Two copy-number states are needed as reference segments; otherwise
    # keep the current solution.
    if (length(lengthlogr) < 2) {
        print("Not possible: fewer than 2 usable copy-number states after applying this shift -- reverting to old solution")
        solution$reverted <- TRUE
        return(solution)
    }
    longest2 <- order(lengthlogr,decreasing=T)[1:2]
    logr1 <- 2^(meanlogr[longest2[1]]/gamma)
    logr2 <- 2^(meanlogr[longest2[2]]/gamma)
    total1 <- as.numeric(names(lengthlogr)[longest2[1]])+shift
    total2 <- as.numeric(names(lengthlogr)[longest2[2]])+shift
    # A reference state can have NA logr (no covered bins); keep the
    # current solution in that case.
    if (is.na(logr1) || is.na(logr2)) {
        print("Not possible: selected reference segment has no underlying data (NA logR) -- reverting to old solution")
        solution$reverted <- TRUE
        return(solution)
    }
    purity <- (2 * logr1/total1/logr2 - 2/total1)/(1 - logr1 *
                                                   total2/logr2/total1 + logr1/logr2/total1 * 2 - 2/total1)
    if(purity>1) purity <- 1
    if(purity<0.1) purity <- 0.1
    ploidy <- (total2 * purity + (1 - purity) * 2)/logr2
    gridpur <- purity + gridpur
    gridpl <- ploidy + gridpl
    # is.finite() also drops NA/NaN/Inf, which a range check alone
    # would keep.
    gridpur <- gridpur[is.finite(gridpur) & gridpur > 0 & gridpur <= 1]
    gridpl <- gridpl[is.finite(gridpl) & gridpl > 0]
    if (length(gridpur) == 0 | length(gridpl) == 0)
    {
        print("Not possible: ploidy<0 or purity \u2209 [0,1] -- reverting to old solution")
        solution$reverted <- TRUE
        return(solution)
    }
    newsol <- searchGrid(track, purs = gridpur, ploidies = gridpl, gamma = gamma,
                         ismale = ismale, isPON = isPON)
    if(newsol$ambiguous)
    {
        print("New solution is ambiguous: reverting to old one")
        solution$reverted <- TRUE
        return(solution)
    }
    newsol$errs <- errs_orig
    newsol
}
