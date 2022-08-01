treatNormal <- function(lCTS, window)
{
    nlCTS <- lapply(lCTS,function(df)
    {
        cat(".")
        nr <- nrow(df)
        starts <- getstartends(end=nr,window=window)$starts
        ends <- getstartends(end=nr,window=window)$ends
        ndf <- data.frame(space=df[starts,"space"],
                          start=df[starts,"start"],
                          end=df[ends,"end"],
                          width=df[starts,"width"]*window,
                          file=df[starts,"file"],
                          records=sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"records"])),
                          nucleotides=sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"nucleotides"])))
        for(i in 8:ncol(df))
            ndf[[colnames(df)[i]]] <- sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],colnames(df)[i]]))
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
