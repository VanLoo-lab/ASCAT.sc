\name{plotSolution}
\alias{plotSolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot fitted profile and underlying data (segmented track,
pre- and post-rounding to integers)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSolution(tracksSingle,purity,ploidy,lwdSeg,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{purity}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{ploidy}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{lwdSeg}{
%%     ~~Describe \code{tracksSingle} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tracksSingle)
{
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
        xlim = c(0, max(breaks)), ylim = c(0, 8), xlab = "Genomic Position",
        ylab = "relative copy number", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    for (i in 1:length(tracksSingle$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
            transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                purity, ploidy), tracksSingle$lCTS[[i]]$end/1e+06 +
                breaks[i], transform_bulk2tumour(tracksSingle$lCTS[[i]]$smoothed,
                purity, ploidy), col = rgb(0.7, 0.7, 0.7, 0.6))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], transform_bulk2tumour(sapply(meansSeg[[i]],
            function(x) x$mu), purity, ploidy), tracksSingle$lSegs[[i]]$output$loc.end/1e+06 +
            breaks[i], transform_bulk2tumour(sapply(meansSeg[[i]],
            function(x) x$mu), purity, ploidy), lwd = lwdSeg,
            col = rgb(0.4, 0.4, 0.4, 0.4))
        segments(tracksSingle$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            tracksSingle$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
            round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            lwd = 2.5, col = rgb(1, 0.5, 0.5))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 2, col = rgb(0.6,
        0.6, 0.6, 0.4))
    text(x = breaks[2:length(breaks)] - 25, y = 5.5, names(breaks)[2:length(breaks)],
        cex = 0.4)
    mtext(side = 3, paste0("guess purity=", purity, "; guess ploidy=",
        ploidy))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
