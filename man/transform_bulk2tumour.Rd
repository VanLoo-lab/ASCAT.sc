\name{transform_bulk2tumour}
\alias{transform_bulk2tumour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform signal from bulk to tumour given purity and ploidy values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_bulk2tumour(x, rho, phi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
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
function (x, rho, phi)
(phi * 2^x - 2 * (1 - rho))/rho
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
