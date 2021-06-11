writeProfile <- function(prof, samplename, outdir)
{
    write.table(prof,
                quote=F,
                sep="\t",
                col.names=T,
                row.names=F,
                file=paste0(outdir,
                            samplename,
                            ".ASCAT.scprofile.txt"))
}
