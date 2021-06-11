\name{fitProfile}
\alias{fitProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit all segmented profile given purity and ploidy
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitProfile(tracksSingle, purity, ploidy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{purity}{
%%     ~~Describe \code{purity} here~~
}
  \item{ploidy}{
%%     ~~Describe \code{ploidy} here~~
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
function (tracksSingle, purity, ploidy)
{
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
                return(list(roundmu = NA, mu = NA, sd = NA, start = out$loc.start[x],
                  end = out$loc.end[x]))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            list(roundmu = transform_bulk2tumour(mu, purity,
                ploidy), mu = mu, sd = sd, start = out$loc.start[x],
                end = out$loc.end[x])
        })
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
