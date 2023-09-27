# ASCAT.sc
Copy number from DNA profiling techniques, including single-cell (sc),
shallow-coverage (sc), and targeted sequencing, as well as methylation arrays


## Install

First download the zip file on disk and then build and install with the following commands

> R CMD BUILD --no-build-vignettes ASCAT.sc/ 

> R CMD INSTALL ASCAT.sc_1.0.tar.gz 
<br>

> R CMD build --no-build-vignettes ASCAT.sc/ 

> R CMD install ASCAT.sc_0.1.tar.gz 
<br>

Alternatively, you can install with devtools in an R session:

> devtools::install_github("VanLoo-lab/ASCAT.sc", build_opts = c("--no-build-vignettes"))

### Dependencies 

The methylation mode now depends on pre-compiled bad loci and panel of
normal data for 450K, Epicv1 (~0.5GB) and Epicv2 arrays,
which are in the R data package *ASCAT.scDataMeth*

As the files are big and github bandwidth is limited, you might want to download the package
(https://drive.google.com/drive/folders/1zDu5-WEYq3OQ8qSZANBOMYTw-H-LNWTU?usp=share_link) as a zip file and
install from R:

> R CMD build ASCAT.scDataMeth

> R CMD install ASCAT.scDataMeth_0.1.tar.gz

Make sure to install dependencies before installing (this might take a while):

> BiocManager::install(c("GenomicRanges", "Biostrings", "DNAcopy", "minfi", "conumee", "Rsamtools", "xgboost"))

#### Epicv2

Install the manifest and annotations:

> BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")

> BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

## Usage

See vignettes for how-to pipelines.


## Cell-line profiles

ASCAT.sc total copy-number profiles on cell lines from various sources (GDSC and CCLE) can be found here:
https://drive.google.com/drive/folders/1RwhF9cw6KP55fHtlNKTzk36syXN9fspo?usp=sharing
These are accompanied by the R script to derive the profiles, as well as sample ID mapping and QC flags.

**Credits to:** Philip S Smith, CRUK Cambridge Institute

