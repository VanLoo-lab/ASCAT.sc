\name{findBestSolution}
\alias{findBestSolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find best solution from the grid
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findBestSolution(tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2, REFs = c(1:23), purs = seq(0.05, 1, 0.005), ploidies = seq(1.7, 6, 0.01), maxTPhi = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{colSeg}{
%%     ~~Describe \code{colSeg} here~~
}
  \item{lwdSeg}{
%%     ~~Describe \code{lwdSeg} here~~
}
  \item{REFs}{
%%     ~~Describe \code{REFs} here~~
}
  \item{purs}{
%%     ~~Describe \code{purs} here~~
}
  \item{ploidies}{
%%     ~~Describe \code{ploidies} here~~
}
  \item{maxTPhi}{
%%     ~~Describe \code{maxTPhi} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2,
    REFs = c(1:23), purs = seq(0.05, 1, 0.005), ploidies = seq(1.7,
        6, 0.01), maxTPhi = 5, ...)
{
    meanS <- mean(unlist(lapply(REFs, function(x) {
        oo <- tracksSingle$lSegs[[x]]$output
        log2(10) * inverse.rle(list(values = oo$seg.mean, lengths = oo$num.mark))
    })), na.rm = T)
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - meanS
        return(x)
    })
    tracksSingle$lSeg <- lapply(tracksSingle$lSeg, function(x) {
        x$output$seg.mean <- log2(10) * x$output$seg.mean - meanS
        return(x)
    })
    searchGrid(tracksSingle, purs, ploidies, maxTumourPhi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
