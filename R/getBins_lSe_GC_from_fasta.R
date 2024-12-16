getBins_lSe_GC_from_fasta <- function(fasta="/PATH/TO/YOUR/REFERENCE/GENOME/fasta.fa",
                                    CHRS = paste0("", c(1:19, "X")), ##contig names
                                    WINDOW=5000,##starting WINDOW = bin size
                                    outdir="./",
                                    filename_save="Build_xx"
                                    )
{
    ## ##########################################
    require(Biostrings)
    ## ##########################################
    REFGENOME <- Biostrings::readDNAStringSet(fasta, format = "fasta")
    print(substr(names(REFGENOME),1,20))
    REFGENOME <- lapply(1:length(CHRS),function(x) REFGENOME[[x]]) #make sure contigs in same order in your fasta as CHRS
    names(REFGENOME) <- CHRS
    allchr <- ALLCHR <- CHRS
    ## ##########################################
    LENGTHCHR <- sapply(REFGENOME,length)
    ## ##########################################
    ## Derive bins
    lSe <- lapply(ALLCHR,function(chr)
    {
        cat(".")
        getStartsEnds(window=WINDOW,
                      chr=chr,
                      lengthChr=LENGTHCHR[chr],
                      dna=REFGENOME)
    })
    names(lSe) <- ALLCHR
    ## ##########################################
    ## Derive GC content within bins
    lGCT <- lapply(ALLCHR,function(chr)
    {
        cat(".")
        gcT <- gcTrack(chr,lSe[[chr]]$starts,
                       lSe[[chr]]$ends,
                       dna=REFGENOME)
    })
    names(lGCT) <- ALLCHR
    ## ##########################################
    save(lGCT,file=paste(outdir,"/lGCT_",filename_save,".",WINDOW,".Rda",sep=""))
    save(lSe,file=paste(outdir,"/lSe_",filename_save,".",WINDOW,".Rda",sep=""))
    return(NULL)
}


## ##########################################
## Test
if(FALSE)
{
    getwd()
    fasta="~/Downloads/GRCm38.primary_assembly.genome.fa.gz"
    CHRS = paste0("", c(1:19, "X")) ##contig names
    WINDOW=5000##starting WINDOW = bin size
    outdir="./"
    filename_save="mm10"
    test=getBins_lSe_GC_from_fasta(fasta,
                                   CHRS=CHRS,
                                   WINDOW=WINDOW,
                                   outdir=outdir,
                                   filename_save=filename_save)
    load(file=paste(outdir,"/lGCT_",filename_save,".",WINDOW,".Rda",sep=""))
    load(file=paste(outdir,"/lSe_",filename_save,".",WINDOW,".Rda",sep=""))
}
## ##########################################


