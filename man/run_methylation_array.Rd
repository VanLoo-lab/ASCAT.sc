\name{run_methylation_array}
\alias{run_methylation_array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to run on methylation arrays (Illumina 450K or EPICv1).
}
\description{
This function runs ASCAT.sc on methylation array data.

Please see vignette and function "run_sc_sequencing" for more info.

}
\usage{
run_methylation_array(idat_dir,
                                  id_normals,
                                  purs = seq(0.1, 1, 0.01),
                                  ploidies = seq(1.7,5, 0.01),
                                  allchr=c(1:22),
                                  maxtumourpsi=5,
                                  segmentation_alpha=0.01,
                                  outdir="./",
                                  conumee=FALSE,
                                  platform=c("450K","Epicv1"),
                                  projectname="project_ma",
                                  predict_refit=TRUE,
                                  print_results=TRUE,
                                  MC.CORES=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id_normals}{
    vector of characters. IDs of normal diploid controls if any in the dataset.
  }
  \item{purs}{
    numeric vector. purity grid values. default suggested for bulk. for
    single cells, try: purs=c(0.5,seq(.99,1,.01)).
}
  \item{ploidies}{
    numeric vector. ploidy grid values. default from 1.7 to 5 can be
    adjusted depending on tumour type. steps of 0.01 can be increased
    for faster runs.
}
  \item{allchr}{
    vector of characters. chromosome names, must be in c(1:22,"X", "Y").
  }
  \item{maxtumourpsi}{
    numeric value. maximum tumour ploidy above which solutions will be
    masked (distance set to infinity in the grid search). helps reduce
    runtime, as some combinations of psi and purs would lead to
    tumour ploidy >> maxtumourpsi.
  }
  \item{segmentation_alpha}{
    segmentation parameter for the circular binary
    segmentation. penalty=1/segmentation_alpha when using multipcf=T
  }
    \item{outdir}{
    character string. path to output directory.
  }
    \item{conumee}{
      boolean. TRUE: use conumee-based binning of probes (if FALSE: do not bin).
  }
  \item{platform}{
    character string. for now, either "450K" or "Epicv1". please contact
    us if you need help with other platforms.
  }
  \item{projectname}{
    character string to use as flag in output file names.
  }
  \item{predict_refit}{
    boolean. use xgboost predictor to predict refitted ploidy and refit the
    profiles automatically.
  }
  \item{print_results}{
    boolean. print results to disk to directory "outdir".
  }
  \item{MC.CORES}{
    numeric value. number of cores to use, relies on parallel's mclapply
  }
}
\details{

}
\value{

  
}
\references{

}
\author{
maxime.tarabichi@ulb.be
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
run_methylation_data; run_target_sequencing;
}
\examples{


## Please see vignette for more examples

####################################################################
idat_dir <- "PATHTODIR/allidat"
id_normals <- c("6042316071_R04C02",
"6042316087_R03C02",
"6042316087_R01C01")
####################################################################


####################################################################
res <- run_methylation_array(idat_dir=idat_dir,
                             id_normals=id_normals,
                             purs = seq(0.05, 1, 0.01),
                             ploidies = seq(1.7,5, 0.01),
                             allchr=paste0("chr",c(1:22)),
                             maxtumourpsi=5,
                             segmentation_alpha=0.01,
                             outdir="testmethyl/",
                             projectname="project_methyl_test",
                             predict_refit=TRUE,
                             print_results=TRUE,
                             MC.CORES=4)
####################################################################


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
