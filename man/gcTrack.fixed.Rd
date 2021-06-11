\name{gcTrack.fixed}
\alias{gcTrack.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derive gc track along a given chromosome (chr) from reference genome
(dna) given bin size (step)

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcTrack.fixed(chr, dna, step = 5e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{dna}{
%%     ~~Describe \code{chr} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
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
function (chr, dna, step = 5e+05)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]], step,
        c("G", "C")))/step
    lGC <- length(gc)
    gc <- gc[seq(step, lGC, step) - round(step/2)]
    names(gc) <- paste("bin", seq(step, lGC, step) - step + 1,
        sep = "-")
    gc
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
