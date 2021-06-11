\name{getTrackForAll}
\alias{getTrackForAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get smoothed coverage track for all chromosomes, starts and ends from bamfile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getTrackForAll(bamfile, window, lSe=NULL, lGCT=NULL, allchr=1:22,
sdNormalise=0, segmentation_alpha=0.01, isDuplicate=F,
isSecondaryAlignment=F, isNotPassingQualityControls=NA,
isUnmappedQuery=NA, mapqFilter=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfile}{
%%     ~~Describe \code{logr} here~~
}
  \item{window}{
%%     ~~Describe \code{CHRS} here~~
}
  \item{lSe}{
%%     ~~Describe \code{STARTS} here~~
}
  \item{lGCT}{
%%     ~~Describe \code{ENDS} here~~
}
  \item{allchr}{
%%     ~~Describe \code{allchr} here~~
}
  \item{sdNormalise}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{segmentation_alpha}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{isDuplicate}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{isSecondaryAlignment}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{isNotPassingQualityControl}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{isUnmappedQuery}{
  %%     ~~Describe \code{allchr} here~~
}
  \item{mapqFilter}{
  %%     ~~Describe \code{allchr} here~~
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
function (logr, CHRS, STARTS, ENDS, allchr = NA)
{
    if (is.na(allchr))
        allchr <- unique(CHRS)
    lT <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        data.frame(start = STARTS[cond], end = ENDS[cond], records = logr[cond],
            fitted = logr[cond], smoothed = logr[cond])
    })
    lSe <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        list(starts = STARTS[cond], ends = ENDS[cond])
    })
    print("segment Tracks")
    lSegs <- lapply(1:length(lT), function(x) {
        require(DNAcopy)
        segments <- segmentTrack(lT[[x]]$smoothed, chr = paste0(x),
            sd = 0, lSe[[x]]$starts, lSe[[x]]$ends)
    })
    names(lSegs) <- paste0(1:length(lT))
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
