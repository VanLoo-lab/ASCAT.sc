\name{sc_plotHeat}
\alias{sc_plotHeat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots all clustered profiles across samples/single-cells as a heatmap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sc_plotHeat(mat, keep1 = rep(T, ncol(mat)), scaleY = 0.3, fundist = "manhattan", funclust = "ward.D2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{keep1}{
%%     ~~Describe \code{keep1} here~~
}
  \item{scaleY}{
%%     ~~Describe \code{scaleY} here~~
}
  \item{fundist}{
%%     ~~Describe \code{fundist} here~~
}
  \item{funclust}{
%%     ~~Describe \code{funclust} here~~
}
  \item{centromeres}{
%%     ~~Describe \code{funclust} here~~
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
function (mat, keep1 = rep(T, ncol(mat)), scaleY = 0.3, fundist = "manhattan",
    funclust = "ward.D2")
{
    require(gplots)
    layout(mat = cbind(c(3, 4), c(2, 1), c(2, 1)), widths = c(0.4,
        1, 1), heights = c(5, 1.4))
    par(mar = c(4, 4, 1, 1))
    xx <- sc_plotGenome(mat, scaleY = scaleY)
    mmm <- t(mat)
    hcc <- hclust(dist(mmm, method = fundist), met = funclust)
    ord <- hcc$order
    mmm[mmm > 1] <- 2
    mmm[mmm < -1] <- -2
    mmm <- mmm[ord[length(ord):1], ]
    im <- array(0.8, dim = c(nrow(mmm), ncol(mmm), 3))
    im[, , 1] <- sc_getCols(mmm, "R")
    im[, , 2] <- sc_getCols(mmm, "G")
    im[, , 3] <- sc_getCols(mmm, "B")
    par(mar = c(0, 4, 0, 1))
    plot(0, 0, col = rgb(0, 0, 0, 0), xlim = c(0, ncol(mmm)),
        ylim = c(0, 1), frame = F, xaxt = "n", yaxt = "n", xlab = "",
        ylab = "")
    rasterImage(im, 0, 0, ncol(mmm), 1)
    abline(v = which(diff(as.numeric(gsub("(.*):(.*)", "\\1",
        colnames(mmm)))) != 0) + 1, col = rgb(1, 1, 1, 0.8))
    axis(side = 2, at = (1:nrow(mmm) - 0.5)/nrow(mmm), rownames(mmm),
        cex.axis = 0.6, col = rgb(0.5, 0.5, 0.5, 0.5), las = 2,
        tick = F, hadj = 0.4)
    par(mar = c(1.5, 0, 1.5, 0))
    gplots:::plot.dendrogram(as.dendrogram(hcc), horiz = TRUE,
        axes = FALSE, yaxs = "i", leaflab = "none")
    plot(0, 0, col = rgb(0, 0, 0, 0), xlim = c(0, 1), ylim = c(0,
        1), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    legend("center", col = c(rgb(sc_getCols(-2, "R"), sc_getCols(-2,
        "G"), sc_getCols(-2, "B")), rgb(sc_getCols(-1, "R"),
        sc_getCols(-1, "G"), sc_getCols(-1, "B")), rgb(sc_getCols(0,
        "R"), sc_getCols(0, "G"), sc_getCols(0, "B")), rgb(sc_getCols(1,
        "R"), sc_getCols(1, "G"), sc_getCols(1, "B")), rgb(sc_getCols(2,
        "R"), sc_getCols(2, "G"), sc_getCols(2, "B"))), legend = c("HD",
        "LOH", "N", "Gain", "Amp"), pch = 19, cex = 1.5, box.col = rgb(0,
        0, 0, 0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
