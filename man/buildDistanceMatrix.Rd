\name{buildDistanceMatrix}
\alias{buildDistanceMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Build distance matrix
}
\description{
Build distance matrix for all combinations of purity and ploidy for the
grid search
}
\usage{
buildDistanceMatrix(meansSeg, weights, purs, ploidies, maxTumourPhi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meansSeg}{
    list object containing mean per segment for each chromosome
}
  \item{weights}{
    weights of all segments (=lengths of the segments in number of data
    points, e.g. bins, probes, usually correlates with genomic lengths)
}
  \item{purs}{
    vector of numeric: realistic purity values to grid (e.g. 0 to 1 by 0.01)
}
  \item{ploidies}{
    vector of numeric: realistic average ploidy values to grid (e.g. 1.1 to 8 by 0.01)
}
  \item{maxTumourPhi}{
    numeric: maximum tumour ploidy to consider (e.g. 8)
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
function (meansSeg, weights, purs, ploidies, maxTumourPhi)
{
    errs <- matrix(NA, length(purs), length(ploidies))
    rownames(errs) <- purs
    colnames(errs) <- ploidies
    for (pp in 1:length(purs)) {
        for (pl in 1:length(ploidies)) {
            if (getTumourPhi(ploidies[pl], purs[pp]) > maxTumourPhi)
                errs[pp, pl] <- Inf
            else errs[pp, pl] <- geterrors(rho = purs[pp], phi = ploidies[pl],
                meansSeg, weights, (pl/pp) * 0 + 1)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
