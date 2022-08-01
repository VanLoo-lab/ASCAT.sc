ts_treatTrack <- function(track,
                       lCTS.normal,
                       lGCT,
                       lSe,
                       lRepli=NULL,
                       window=ceiling(as.numeric(WINDOW)/as.numeric("10000")))
{
    lCTS <- lapply(track$lCTS,function(df)
    {
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
    })
    print("correct for GC, replication timing and other biases in normal")
    lCTS <- ts_smoothCoverageTrackAll(lCTS,
                                      lSe = lSe,
                                      lGCT = lGCT,
                                      lRepli = lRepli,
                                      lNormal = lCTS.normal)
    print("segment Tracks")
    lSegs <- lapply(1:length(lCTS), function(x) {
        require(DNAcopy)
        segments <- segmentTrack(lCTS[[x]]$smoothed,
                                 chr = paste0(x),
                                 sd = 0,
                                 lSe[[x]]$starts,
                                 lSe[[x]]$ends)
    })
    names(lSegs) <- paste0(1:length(lCTS))
    tracks <- list(lCTS = lCTS, lSegs = lSegs)
    return(tracks)
}
