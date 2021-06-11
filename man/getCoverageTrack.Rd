\name{getCoverageTrack}
\alias{getCoverageTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get coverage track from BAM file (bamPath) along a given chromosome
(chr) and starts/ends of bins
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCoverageTrack(bamPath, chr, starts, ends, CHRSTRING = "", isDuplicate = F, isSecondaryAlignment = F, isNotPassingQualityControls = NA, isUnmappedQuery = NA, mapqFilter = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamPath}{
%%     ~~Describe \code{bamPath} here~~
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
  \item{CHRSTRING}{
%%     ~~Describe \code{CHRSTRING} here~~
}
  \item{isDuplicate}{
%%     ~~Describe \code{isDuplicate} here~~
}
  \item{isSecondaryAlignment}{
%%     ~~Describe \code{isSecondaryAlignment} here~~
}
  \item{isNotPassingQualityControls}{
%%     ~~Describe \code{isNotPassingQualityControls} here~~
}
  \item{isUnmappedQuery}{
%%     ~~Describe \code{isUnmappedQuery} here~~
}
  \item{mapqFilter}{
%%     ~~Describe \code{mapqFilter} here~~
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
function (bamPath, chr, starts, ends, CHRSTRING = "", isDuplicate = F,
    isSecondaryAlignment = F, isNotPassingQualityControls = NA,
    isUnmappedQuery = NA, mapqFilter = 0)
{
    require(Rsamtools)
    require(GenomicRanges)
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = isDuplicate,
        isSecondaryAlignment = isSecondaryAlignment, isNotPassingQualityControls = isNotPassingQualityControls,
        isUnmappedQuery = isUnmappedQuery), which = GRanges(paste0(CHRSTRING,
        chr), IRanges(starts, ends)), mapqFilter = mapqFilter)
    coverageTrack <- countBam(bamPath, param = sbp)
    coverageTrack$records <- coverageTrack$records
    return(coverageTrack)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
