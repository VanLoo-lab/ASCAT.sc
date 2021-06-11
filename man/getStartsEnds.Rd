\name{getStartsEnds}
\alias{getStartsEnds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get pre-computed starts and ends for a given window and chromosome
extracted from disk (pathWindows) or usig fixed steps from 0 to
lengthChr by window. Optionally can remove badbins from pre-computed track.
}
\description{
Get pre-computed starts and ends for a given window and chromosome
extracted from disk (pathWindows) or usig fixed steps from 0 to
lengthChr by window. Optionally can remove badbins from pre-computed track.
}
\usage{
getStartsEnds(window, chr, pathWindows = NA, pathBadBins = NA, lengthChr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{window}{
%%     ~~Describe \code{window} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{pathWindows}{
%%     ~~Describe \code{pathWindows} here~~
}
  \item{pathBadBins}{
%%     ~~Describe \code{pathBadBins} here~~
}
  \item{lengthChr}{
%%     ~~Describe \code{lengthChr} here~~
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
function (window, chr, pathWindows = NA, pathBadBins = NA, lengthChr)
{
    if (!is.na(pathWindows)) {
        badbins <- if (!is.na(pathBadBins))
            read.table(pathBadBins)[, 1]
        else NULL
        t <- read.table(pathWindows, header = T)
        t <- t[-c(badbins), ]
        subt <- t[t$CHR == chr, ]
        starts <- c(1, as.numeric(as.character(subt$END[-c(length(subt$END))])) +
            1)
        ends <- as.numeric(as.character(subt$END))
    }
    else {
        divideChr <- seq(0, lengthChr, window)
        starts <- divideChr[-c(length(divideChr))] + 1
        ends <- divideChr[-c(1)]
    }
    list(starts = starts, ends = ends)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
