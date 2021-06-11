getStartsEnds <- function(window,
                          chr,
                          lengthChr,
                          dna=NULL,
                          pathWindows=NA,
                          pathBadBins=NA,
                          centromeres=NULL,
                          excludebed=NULL)
{
    if(!is.na(pathWindows))
    {
        badbins <- if(!is.na(pathBadBins)) read.table(pathBadBins)[,1] else NULL
        t <- read.table(pathWindows,header=T)
        t <- t[-c(badbins),]
        subt <- t[t$CHR==chr,]
        starts <- c(1,as.numeric(as.character(subt$END[-c(length(subt$END))]))+1)
        ends <- as.numeric(as.character(subt$END))
    }
    else if(!is.null(excludebed))
    {
        keepChr <- gsub("chr","",excludebed[,1])==gsub("chr","",chr)
        if(sum(keepChr)>0)
            {
                excludebed <- excludebed[keepChr,,drop=F]
                excludebed <- rbind(excludebed,excludebed[nrow(excludebed),])
                excludebed[nrow(excludebed),2:3] <- lengthChr-1:0
                grE <- GRanges(chr,IRanges(excludebed[,2],excludebed[,3]))
                covs <- coverage(grE)
                effectiveLength <- sum(covs==0)
                divideChr <- seq(1, effectiveLength, window)
                starts <- divideChr[-c(length(divideChr))]
                effectivePos <- which(as.logical((covs==0)[[1]]))
                starts <- effectivePos[starts]
                ends <- c(starts[-c(1)]-1,lengthChr)
            }
        else
        {
            return(getStartsEnds(window=window,
                                 chr=chr,
                                 lengthChr=lengthChr,
                                 centromeres=centromeres))
        }
    }
    else
    {
        divideChr <- seq(0, lengthChr, window)
        starts <- divideChr[-c(length(divideChr))] + 1
        ends <- divideChr[-c(1)]
    }
    excludeBadBins(removeCentromeres(list(starts=starts,ends=ends),chr=chr,centromeres=centromeres),
                   chr=chr,dna=dna)
}
