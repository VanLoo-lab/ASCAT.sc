\name{smoothCoverageTrack}
\alias{smoothCoverageTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
smooth coverage track
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smoothCoverageTrack(lCT, lSe, lGCT, method = c("loess", "lowess"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lCT}{
%%     ~~Describe \code{lCT} here~~
}
  \item{lSe}{
%%     ~~Describe \code{lSe} here~~
}
  \item{lGCT}{
%%     ~~Describe \code{lGCT} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (lCT, lSe, lGCT, method = c("loess", "lowess"))
{
    allRec <- unlist(lapply(lCT, function(x) log10(x$records +
        1)))
    allGC <- unlist(lGCT)
    starts <- c(0, cumsum(sapply(lCT, nrow)[-c(length(lCT))])) +
        1
    ends <- cumsum(sapply(lCT, nrow))
    if (method[1] == "lowess")
        smoothT <- lowess(allGC, allRec, f = 2/3, iter = 10)
    if (method[1] == "loess")
        smoothT <- loess(allRec ~ allGC)
    for (i in 1:length(lCT)) {
        lCT[[i]] <- cbind(lCT[[i]], smoothT$fitted[starts[i]:ends[i]],
            smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]]) - 1):ncol(lCT[[i]])] <- c("fitted",
            "smoothed")
    }
    return(lCT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
