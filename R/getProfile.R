getProfile <- function(profile, CHRS=c(1:22))
{
    names(profile) <- CHRS
    tt <- NULL
    for(i in 1:length(profile))
    {
        for(j in 1:length(profile[[i]]))
        {
            tt <- rbind(tt,c(CHRS[i],
                        profile[[i]][[j]]$start,
                        profile[[i]][[j]]$end,
                        profile[[i]][[j]]$num.mark,
                        round(profile[[i]][[j]]$roundmu),
                        (profile[[i]][[j]]$roundmu),
                        (profile[[i]][[j]]$mu),
                        (profile[[i]][[j]]$sd)))
        }
    }
    colnames(tt) <- c("chromosome",
                      "start",
                      "end",
                      "num.mark",
                      "total_copy_number",
                      "total_copy_number_logr",
                      "logr",
                      "logr.sd")
    tt <- data.frame(chromosome=as.character(tt[,"chromosome"]),
                     start=as.numeric(tt[,"start"]),
                     end=as.numeric(tt[,"end"]),
                     num.mark=as.numeric(tt[,"num.mark"]),
                     total_copy_number=as.numeric(tt[,"total_copy_number"]),
                     total_copy_number_logr=as.numeric(tt[,"total_copy_number_logr"]),
                     logr=as.numeric(tt[,"logr"]),
                     logr.sd=as.numeric(tt[,"logr.sd"]))
    tt
}
