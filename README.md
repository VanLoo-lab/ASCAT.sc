# ASCAT.sc
Copy number from DNA profiling techniques, including single-cell (sc),
shallow-coverage (sc), and targeted sequencing, as well as methylation arrays


## Install

You can install from github with the following command:

```{r}
devtools::install_github("VanLoo-lab/ASCAT.sc", build_opts = c("--no-build-vignettes"))
```

### Dependencies 

Make sure to install dependencies before installing (this might take a while):

```{r}
devtools::install_github("iovlaicu/copynumber")
BiocManager::install(c("GenomicRanges", "Biostrings", "DNAcopy", "Rsamtools", "xgboost"))
```


## Examples and Usage

Go to [Wiki](https://github.com/VanLoo-lab/ASCAT.sc/wiki) pages.

Vignettes are now obsolete, but remain a good way to understand the
different steps of the pipelines.

## Cell-line profiles

ASCAT.sc total copy-number profiles on cell lines from various sources (GDSC and CCLE) can be found [here](https://drive.google.com/drive/folders/1RwhF9cw6KP55fHtlNKTzk36syXN9fspo?usp=sharing).
These are accompanied by the R script to derive the profiles, as well as sample ID mapping and QC flags.

**Credits to:** Philip S Smith, CRUK Cambridge Institute


## See also

ASCAT.ma is a separate package, that is built in ASCAT.sc' framework
to derive genomewide copy-number calls, purity and ploidy estimates
from (Illumina) methylation arrays.
