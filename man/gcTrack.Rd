\name{gcTrack}
\alias{gcTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derive gc track along a given chromosome (chr) from reference genome
(dna) given starts and ends of bins
}
\description{
Derive gc track along a given chromosome (chr) from reference genome
(dna) given starts and ends of bins
}
\usage{
gcTrack(chr, starts, ends, dna, window = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{ends}{
%%     ~~Describe \code{ends} here~~
}
  \item{dna}{
%%     ~~Describe \code{dna} here~~
}
  \item{window}{
    smooth window for moving average of GC content along genome
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
function (chr, starts, ends, dna, window = 5000)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]], window,
        c("G", "C")))/window
    if (ends[length(ends)] > length(gc)) {
        ends[length(ends)] <- length(gc)
    }
    gc <- sapply(1:length(starts), function(x) mean(gc[starts[x]:ends[x]]))
    names(gc) <- paste("bin", 1:length(starts), sep = "-")
    gc
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
