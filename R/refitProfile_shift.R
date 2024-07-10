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
    keep <- rep(T,length(lengthlogr))
    if(shift==-1) {keep <- !names(lengthlogr)%in%c("0","1")}
    lengthlogr <- lengthlogr[keep]
    meanlogr <- meanlogr[keep]
    longest2 <- order(lengthlogr,decreasing=T)[1:2]
    logr1 <- 2^(meanlogr[longest2[1]]/gamma)
    logr2 <- 2^(meanlogr[longest2[2]]/gamma)
    total1 <- as.numeric(names(lengthlogr)[longest2[1]])+shift
    total2 <- as.numeric(names(lengthlogr)[longest2[2]])+shift
    purity <- (2 * logr1/total1/logr2 - 2/total1)/(1 - logr1 *
                                                   total2/logr2/total1 + logr1/logr2/total1 * 2 - 2/total1)
    if(purity>1) purity <- 1
    if(purity<0.1) purity <- 0.1
    ploidy <- (total2 * purity + (1 - purity) * 2)/logr2
    gridpur <- purity + gridpur
    gridpl <- ploidy + gridpl
    gridpur <- gridpur[gridpur > 0 & gridpur <= 1]
    gridpl <- gridpl[gridpl > 0]
    if (length(gridpur) == 0 | length(gridpl) == 0)
    {
        print("Not possible: ploidy<0 or purity ∉ [0,1]")
        return(list(errs=NULL,
                    purity=NA,
                    ploidy=NA,
                    ambiguous=TRUE))
    }
    newsol <- searchGrid(track, purs = gridpur, ploidies = gridpl, gamma = gamma,
                         ismale = ismale, isPON = isPON)
    if(newsol$ambiguous)
    {
        print("New solution is ambiguous: reverting to old one")
        solution$reverted <- TRUE
        return(solution)
    }
    newsol
}
##########################################################################
