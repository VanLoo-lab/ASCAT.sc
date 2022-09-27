myloess <- function(LL, ...)
{
    if(LL>10000)
        stats::loess(...,
                     control=loess.control(surface = "interpolate",
                                           statistics = "approximate",
                                           trace.hat = "approximate",
                                           cell = 0.2,
                                           iterations = 4,
                                           iterTrace = FALSE))
    else
        stats::loess(...)
}
