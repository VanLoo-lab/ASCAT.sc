sc_plotHeat <- function(mat,
                        keep1=rep(T,ncol(mat)),
                        scaleY=.3,
                        fundist="manhattan",
                        funclust="ward.D2",
                        centromeres=NULL)
{
    require(gplots)
    require(rdist)
    layout(mat=cbind(c(3,4),c(2,1),c(2,1)),
           widths=c(.4,1,1),
           heights=c(5,1.4))
    par(mar=c(4,4,1,1))
    xx <- sc_plotGenome(mat,scaleY=scaleY,centromeres=centromeres)
    mmm <- t(mat)
    hcc <- hclust(rdist(mmm,metric=fundist),met=funclust)
    ord <- hcc$order
    mmm[mmm>1] <- 2
    mmm[mmm< -1] <- -2
    mmm <- mmm[ord[length(ord):1],]
    im <- array(0.8,dim=c(nrow(mmm),ncol(mmm),3))
    im[,,1] <- sc_getCols(mmm,"R")
    im[,,2] <- sc_getCols(mmm,"G")
    im[,,3] <- sc_getCols(mmm,"B")
    par(mar=c(0,4,0,1))
    plot(0,0,col=rgb(0,0,0,0),xlim=c(0,ncol(mmm)),ylim=c(0,1),frame=F,
         xaxt="n",yaxt="n",xlab="",ylab="")
    rasterImage(im,0,0,ncol(mmm),1)
    abline(v=which(diff(as.numeric(gsub("(.*):(.*)","\\1",
                                        colnames(mmm))))!=0)+1,
           col=rgb(1,1,1,.8))
    axis(side=2,at=(1:nrow(mmm)-.5)/nrow(mmm),rownames(mmm),
         cex.axis=.6,col=rgb(.5,.5,.5,.5),las=2,tick=F,hadj=.4)
    par(mar=c(1.5,0,1.5,0))
    gplots:::plot.dendrogram(as.dendrogram(hcc),
                             horiz = TRUE,
                             axes = FALSE,
                             yaxs = "i",
                             edgePar=list(t.cex=1/log10(ncol(mat))),
                             leaflab = "none")
    plot(0,0,col=rgb(0,0,0,0),
         xlim=c(0,1),
         ylim=c(0,1),
         frame=F,
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="")
    legend("center",col=c(rgb(sc_getCols(-2,"R"),sc_getCols(-2,"G"),
                              sc_getCols(-2,"B")),
                          rgb(sc_getCols(-1,"R"),sc_getCols(-1,"G"),
                              sc_getCols(-1,"B")),
                          rgb(sc_getCols(0,"R"),sc_getCols(0,"G"),sc_getCols(0,"B")),
                          rgb(sc_getCols(1,"R"),sc_getCols(1,"G"),sc_getCols(1,"B")),
                          rgb(sc_getCols(2,"R"),sc_getCols(2,"G"),sc_getCols(2,"B"))),
           legend=c("HD","LOH","N","Gain","Amp"),
           pch=19,cex=1.5,box.col=rgb(0,0,0,0))
}
