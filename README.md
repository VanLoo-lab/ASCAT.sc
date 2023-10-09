# ASCAT.sc
Copy number from DNA profiling techniques, including single-cell (sc),
shallow-coverage (sc), and targeted sequencing, as well as methylation arrays


## Install

First download the zip file on disk and then build and install with the following commands

```{r}
devtools::install_github("VanLoo-lab/ASCAT.sc", build_opts = c("--no-build-vignettes"))
```

### Dependencies 

Make sure to install dependencies before installing (this might take a while):

```{r}
BiocManager::install(c("GenomicRanges", "Biostrings", "DNAcopy", "minfi", "conumee", "Rsamtools", "xgboost"))
```

#### Methylation arrays
The methylation mode now depends on pre-compiled bad probes, and panel of
normal diploid data for 450K, Epicv1 (~0.5GB) and Epicv2 arrays,
which are in the R data package *ASCAT.scDataMeth*.

As the files are big and github bandwidth is limited, you can download
the package as a zip file
[here]((https://drive.google.com/drive/folders/1zDu5-WEYq3OQ8qSZANBOMYTw-H-LNWTU?usp=share_link))
and build/install for R as follows:

> R CMD build ASCAT.scDataMeth

> R CMD install ASCAT.scDataMeth_0.1.tar.gz

For Epicv2, also install the manifest and annotations:

```{r}
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
```

## Examples and Usage

Go to [Wiki](https://github.com/VanLoo-lab/ASCAT.sc/wiki) pages.

Vignettes are now obsolete, but remain a good way to understand the
different steps of the pipelines.

## Cell-line profiles

ASCAT.sc total copy-number profiles on cell lines from various sources (GDSC and CCLE) can be found [here](https://drive.google.com/drive/folders/1RwhF9cw6KP55fHtlNKTzk36syXN9fspo?usp=sharing).
These are accompanied by the R script to derive the profiles, as well as sample ID mapping and QC flags.

**Credits to:** Philip S Smith, CRUK Cambridge Institute

