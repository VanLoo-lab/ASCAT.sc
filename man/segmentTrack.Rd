\name{segmentTrack}
\alias{segmentTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segments coverage track
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segmentTrack(covtrack, chr, starts, ends = NA, sd = 0, min.width = 5, ALPHA = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covtrack}{
%%     ~~Describe \code{covtrack} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{ends}{
%%     ~~Describe \code{ends} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{min.width}{
%%     ~~Describe \code{min.width} here~~
}
  \item{ALPHA}{
%%     ~~Describe \code{ALPHA} here~~
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
function (covtrack, chr, starts, ends = NA, sd = 0, min.width = 5,
    ALPHA = 0.01)
{
    covtrack <- covtrack * rnorm(length(covtrack), mean = 1,
        sd = sd)
    cna <- CNA(covtrack, chr = rep(chr, length(covtrack)), maploc = starts,
        data.type = "logratio")
    cna.smoothed <- smooth.CNA(cna)
    segment(cna.smoothed, min.width = min.width, alpha = ALPHA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
