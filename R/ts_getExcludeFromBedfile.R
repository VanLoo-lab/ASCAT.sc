ts_getExcludeFromBedfile <- function(bed,chr)
{
    subset <- gsub("chr","",bed[,1])==chr
    if(sum(subset)>0)
        return(list(starts=bed[subset,2],
                    ends=bed[subset,3]))
    else return(list(starts=1,
                     ends=2))
}
