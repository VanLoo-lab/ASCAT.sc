\name{writeProfile}
\alias{writeProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Writes profile on disk
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeProfile(prof, samplename, outdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prof}{
%%     ~~Describe \code{prof} here~~
}
  \item{samplename}{
%%     ~~Describe \code{samplename} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
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
function (prof, samplename, outdir)
{
    write.table(prof, quote = F, sep = "\t", col.names = T, row.names = F,
        file = paste0(outdir, samplename, ".maxPhiT5.profile.txt"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
