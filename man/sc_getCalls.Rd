\name{sc_getCalls}
\alias{sc_getCalls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gathers all calls from segmented and fitted profiles of single cells or
related samples.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sc_getCalls(tracksSingle, mc.cores = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
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
function (tracksSingle, mc.cores = 10)
{
    require(parallel)
    medianS <- median(unlist(lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed
    })))
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - medianS
        return(x)
    })
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- mclapply(1:nrow(out), function(x) {
            cond1 <- tracksSingle$lCTS[[i]]$start >= out$loc.start[x]
            cond2 <- tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            isIn <- cond1 & cond2
            if (sum(isIn) < 2)
                return(list(roundmu = NA, mu = NA, sd = NA))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            list(roundmu = round(2^(mu + log2(2)) - 2), mu = mu,
                sd = sd)
        }, mc.cores = mc.cores)
        t <- cbind(out, sapply(means, function(y) y$roundmu),
            sapply(means, function(y) y$mu), sapply(means, function(y) y$sd))
        colnames(t) <- c(colnames(out), "CNtot", "CNseg", "CNsd")
        t
    })
    return(meansSeg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
