ts_smoothTrack <- function(R,GC, NORMAL=NULL, REPLI=NULL)
{
    if(!is.null(NORMAL))
        R <- R-NORMAL
    ## credit @Jonas Demeulemeester for the spline correction
    corrdata <- data.frame(r = R,
                           GC = GC,
                           repli = if(is.null(REPLI)) rep(NA,length(GC)) else REPLI)
    model <- if(is.null(REPLI))
                     lm(r ~ splines::ns(x = GC, df = 5, intercept = T),
                        y=F,
                        model = F,
                        data = corrdata,
                        na.action="na.exclude")
             else
                 lm(r ~ splines::ns(x = GC, df = 5, intercept = T)
                    + splines::ns(x = repli, df = 5, intercept = T),
                    y=F,
                    model = F,
                    data = corrdata,
                    na.action="na.exclude")
    list(residuals=residuals(model),fitted=fitted(model))
}
