\name{isGoodSolution}
\alias{isGoodSolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks if a fitted solution is acceptable (i.e. the size of homozygous
deletion is below a given threshold in megabases).
}
\description{
Checks if a fitted solution is acceptable (i.e. the size of homozygous
deletion is below a given threshold in megabases).
}
\usage{
isGoodSolution(meansSeg, maxSize0 = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meansSeg}{
%%     ~~Describe \code{meansSeg} here~~
}
  \item{maxSize0}{
%%     ~~Describe \code{maxSize0} here~~
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
function (meansSeg, maxSize0 = 20)
{
    starts <- unlist(lapply(meansSeg, function(x) lapply(x, function(y) y$start)))
    ends <- unlist(lapply(meansSeg, function(x) lapply(x, function(y) y$end)))
    sizes <- (ends - starts)/1e+06
    cn <- unlist(lapply(meansSeg, function(x) lapply(x, function(y) y$roundmu)))
    cn0 <- round(cn) <= 0
    totalsizes0 <- sum(sizes[cn0], na.rm = T)
    if (totalsizes0 < maxSize0)
        return(TRUE)
    FALSE
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
