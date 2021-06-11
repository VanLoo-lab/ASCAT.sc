\name{sc_getCols}
\alias{sc_getCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns colour intensity for a given channel and a given copy number (mmm)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sc_getCols(mmm, channel = "R")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mmm}{
%%     ~~Describe \code{mmm} here~~
}
  \item{channel}{
%%     ~~Describe \code{channel} here~~
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
function (mmm, channel = "R")
{
    poss <- as.character(c(-2, -1, 0, 1, 2))
    Rs <- c(0.6, 0.3, 0, 0, 0)
    Gs <- c(0, 0, 0, 0.3, 0.6)
    Bs <- c(0.2, 0.2, 0, 0.5, 0.5)
    names(Rs) <- names(Gs) <- names(Bs) <- poss
    get(paste0(channel, "s"))[as.character(mmm)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
