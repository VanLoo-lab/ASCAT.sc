ts_removeOnTargets <- function(lCT,
                               lCTex)
{
    for(i in 1:length(lCT))
    {
        gr1 <- GRanges(lCT[[i]][,1],IRanges(lCT[[i]][,2],lCT[[i]][,3]))
        gr2 <- GRanges(lCTex[[i]][,1],IRanges(lCTex[[i]][,2],lCTex[[i]][,3]))
        ov <- findOverlaps(gr2,gr1, type="within")
        rem <- tapply(1:length(queryHits(ov)),subjectHits(ov),function(x) sum(lCTex[[i]][queryHits(ov)[x],"records"],na.rm=T))
        lCT[[i]][as.numeric(names(rem)),"records"] <- lCT[[i]][as.numeric(names(rem)),"records"]-rem
    }
    lCT
}
