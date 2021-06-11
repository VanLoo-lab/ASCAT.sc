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
                        round(profile[[i]][[j]]$roundmu),
                        (profile[[i]][[j]]$roundmu),
                        (profile[[i]][[j]]$mu),
                        (profile[[i]][[j]]$sd)))
        }
    }
    colnames(tt) <- c("chromosome",
                      "start",
                      "end",
                      "total_copy_number",
                      "total_copy_number_logr",
                      "logr",
                      "logr.sd")
    tt
}
