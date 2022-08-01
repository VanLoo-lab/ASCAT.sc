getstartends <- function(start=1, end, window)
{
    starts <- seq(start,end,window)
    starts <- starts[-c(length(starts))]
    ends <- c(starts[-c(1)]-1,end)
    list(starts=starts,ends=ends)
}
