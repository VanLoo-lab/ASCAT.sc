\name{getProfile}
\alias{getProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get absolute copy number number profile from segmented smoothed and
fitted track (profile)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getProfile(profile, CHRS = c(1:22))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{profile}{
    fitted profile
}
  \item{CHRS}{
%%     ~~Describe \code{CHRS} here~~
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
function (profile, CHRS = c(1:22))
{
    names(profile) <- CHRS
    tt <- NULL
    for (i in 1:length(profile)) {
        for (j in 1:length(profile[[i]])) {
            tt <- rbind(tt, c(CHRS[i], profile[[i]][[j]]$start,
                profile[[i]][[j]]$end, round(profile[[i]][[j]]$roundmu),
                (profile[[i]][[j]]$roundmu), (profile[[i]][[j]]$mu),
                (profile[[i]][[j]]$sd)))
        }
    }
    colnames(tt) <- c("chromosome", "start", "end", "total_copy_number",
        "total_copy_number_logr", "logr", "logr.sd")
    tt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
