\name{run_sc_sequencing}
\alias{run_sc_sequencing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to run on shallow-coverage/single-cell (sc) sequencing data.
}
\description{
This function bins the genome into small bins and keeps the ones with
high consistency of depth of coverage in diploid samples. It counts
(high mapping quality and properly mapped) reads/pairs of reads within
those bins then applies the ASCAT methodology for smoothing (GC
correction with loess), segmenting (CBS or multipcf for multi-region sampling) and fitting of purity/ploidy and integer
copy-number.
It will also automatically refit the profiles to match the feature space
of high-quality consensus PCAWG calls.
Finally it prints the results and profiles to disk and returns a list with all raw and
processed data.
  
}
\usage{
run_sc_sequencing(tumour_bams,
                              allchr=paste0("",c(1:22,"X")),
                              sex=c("female","male","female"),
                              chrstring_bam="",
                              purs = seq(0.1, 1, 0.01),
                              ploidies = seq(1.7,5, 0.01),
                              maxtumourpsi=5,
                              binsize=500000,
                              segmentation_alpha=0.01,
                              predict_refit=TRUE,
                              print_results=TRUE,
                              build=c("hg19","hg38"),
                              MC.CORES=1,
                              outdir="./",
                              probs_filters=.1,
                              sc_filters=FALSE,
                              projectname="project",
                              multipcf=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumour_bams}{
    vector of characters. paths to the bams (sorted and indexed).
}
  \item{allchr}{
    vector of characters. chromosome names, must be in c(1:22,"X", "Y"),
    can have a "chr" prefix.
}
  \item{sex}{
    vector of characters. sex in "male" or "female" for each of the
    patients with tumours in tumour_bams.
}
  \item{chrstring_bam}{
    character. either "" or "chr" corresponding to presence of string in
    contig names in the bam.
}
  \item{purs}{
    numeric vector. purity grid values.
}
  \item{ploidies}{
    numeric vector. ploidy grid values.
}
  \item{maxtumourpsi}{
    numeric value. maximum tumour ploidy above which solutions will be
    masked (distance set to infinity in the grid search).
}
  \item{binsize}{
    the size of the bins, should be above 30kb and will be set to a multiple of 30kb, which is
    the current smallest bin size.
}
  \item{segmentation_alpha}{
    segmentation parameter for the circular binary
    segmentation. penalty=1/segmentation_alpha when using multipcf=T
}
  \item{predict_refit}{
    boolean. use xgboost predictor to predict refitted ploidy and refit the
    profiles automatically.
}
  \item{print_results}{
    boolean. print results to disk.
}
  \item{build}{
    character string. either "hg19" or "hg38".
}
  \item{MC.CORES}{
    numeric value. number of cores to use, relies on parallel's mclapply
}
  \item{outdir}{
    character string. path to output directory.
  }
    \item{probs_filters}{
      numeric value. if sc_filters is TRUE, proportion of bad cells
      expected. leave at .1, more will be identified and removed anyways.
  }
    \item{sc_filters}{
      boolean. for single-cell data, should bad cells be identified and
      flagged. will generate a pdf with quality checks.
}
  \item{projectname}{
    character string to use as flag in output file names.
}
  \item{multipcf}{
    boolean. use multi track segmentation when bams come from regions of
    the same tumour (i.e. breakpoints are likely shared).
}
}
\details{

}
\value{

  
}
\references{

}
\author{

}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
run_methylation_data; run_target_sequencing;
}
\examples{




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
