getRefGenome <-
function(fasta=FASTA,
         CHRS=paste0("",c(1:22,"X","Y","MT")))
{
    dna <- Biostrings::readDNAStringSet(fasta,format="fasta")
    dna <- lapply(1:length(CHRS),function(x) dna[[x]])
    names(dna) <- CHRS
    return(dna)
}
