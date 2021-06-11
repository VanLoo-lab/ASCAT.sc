\name{sc_plotGenome}
\alias{sc_plotGenome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot copy number profiles of series of samples/single cells along the genome
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sc_plotGenome(mat, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2, scaleY = 0.6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{colSeg}{
%%     ~~Describe \code{colSeg} here~~
}
  \item{lwdSeg}{
%%     ~~Describe \code{lwdSeg} here~~
}
  \item{scaleY}{
%%     ~~Describe \code{scaleY} here~~
}
  \item{centromeres}{
%%     ~~Describe \code{scaleY} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (mat, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2,
    scaleY = 0.6, ...)
{
    centromeres[, 1] <- gsub("chr", "", centromeres[, 1])
    matPos <- mat
    matPos[matPos < 0] <- 0
    matPos[matPos > 0] <- 1
    matNeg <- mat
    matNeg[matNeg > 0] <- 0
    matNeg[matNeg < 0] <- -1
    vecPositive <- rowSums(matPos)
    names(vecPositive) <- rownames(mat)
    vecNegative <- rowSums(matNeg)
    names(vecNegative) <- rownames(mat)
    par(mar = c(4, 4, 1, 1))
    chrom <- as.numeric(gsub("(.*):(.*)-(.*)", "\\\\1", names(vecPositive)))
    starts <- as.numeric(gsub("(.*):(.*)-(.*)", "\\\\2", names(vecPositive)))
    ends <- as.numeric(gsub("(.*):(.*)-(.*)", "\\\\3", names(vecPositive)))
    breaks <- c(0, cumsum(sapply(unique(chrom), function(x) max(ends[chrom ==
        x]))))/1e+06
    names(breaks) <- c(0, unique(chrom))
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
        xlim = c(0, max(breaks)), ylim = c(-ncol(mat), ncol(mat)) *
            scaleY, xlab = "Cumulative Genomic Position (Mb)",
        ylab = "Number of Samples", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    for (i in 1:length(unique(chrom))) {
        takeChrom <- chrom == unique(chrom)[i]
        segments(starts[takeChrom]/1e+06 + breaks[i], vecPositive[takeChrom],
            ends[takeChrom]/1e+06 + breaks[i], vecPositive[takeChrom],
            col = rgb(getCols(1, "R"), getCols(1, "G"), getCols(1,
                "B")), lwd = lwdSeg)
        segments(starts[takeChrom]/1e+06 + breaks[i], vecNegative[takeChrom],
            ends[takeChrom]/1e+06 + breaks[i], vecNegative[takeChrom],
            col = rgb(getCols(-1, "R"), getCols(-1, "G"), getCols(-1,
                "B")), lwd = lwdSeg)
        abline(v = centromeres[centromeres[, 1] == i, 2]/1e+06 +
            breaks[i], lty = 2, col = rgb(0.4, 0.4, 0.4, 0.4))
        abline(v = centromeres[centromeres[, 1] == i, 3]/1e+06 +
            breaks[i], lty = 2, col = rgb(0.4, 0.4, 0.4, 0.4))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 1, col = rgb(0.2,
        0.1, 0.1, 1))
    text(x = breaks[2:length(breaks)] - 25, y = ncol(mat) * scaleY,
        names(breaks)[2:length(breaks)], cex = 0.4)
    max(breaks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
