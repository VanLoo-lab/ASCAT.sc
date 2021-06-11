\name{getCoverageTrack.Fix}
\alias{getCoverageTrack.Fix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get coverage track from BAM file (bamPath) along a given chromosome
(chr), the length of the chromosome (lengthChr) and the fixed bin size (step)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCoverageTrack.Fix(bamPath, chr, lengthChr, step, CHRSTRING = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamPath}{
%%     ~~Describe \code{bamPath} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{lengthChr}{
%%     ~~Describe \code{lengthChr} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
}
  \item{CHRSTRING}{
%%     ~~Describe \code{CHRSTRING} here~~
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
function (bamPath, chr, lengthChr, step, CHRSTRING = "")
{
    require(Rsamtools)
    require(GenomicRanges)
    divideChr <- seq(0, lengthChr, step)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE),
        which = GRanges(paste0(CHRSTRING, chr), IRanges(starts,
            ends)))
    coverageTrack <- countBam(bamPath, param = sbp)
    return(coverageTrack)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
