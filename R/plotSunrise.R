plotSunrise <- function(solution)
{
    plot(0,0,col=rgb(0,0,0,0),
         xlab="ploidy",
         ylab="purity",
         xaxt="n",yaxt="n",frame=F,xlim=c(0,1),ylim=c(0,1))
    errs <- solution$errs
    errs <- errs-min(errs)
    errs.max <- max(solution$errs[!is.infinite(solution$errs)])
    errs[is.infinite(errs)] <- errs.max
    errs <- errs/errs.max
    .getCol <- function(x)
    {
        vec <- c(sapply(seq(0,1,length.out=100),function(x) rgb(.8,0,0,1-x)),
                 sapply(seq(0,1,length.out=100),function(x) rgb(0,0,.8,x)))
        vec[pmax(1,pmin(200,round(x*200)))]
    }
    im <- matrix(.getCol(as.vector(1-errs)),dim(errs)[1],dim(errs)[2])
    rasterImage(as.raster(im),0,0,1,1)
    points(which(colnames(errs)==solution$ploidy)/ncol(errs),
           1-which(rownames(errs)==solution$purity)/nrow(errs),
           col="green",pch=19)
    axis(side=1,at=seq(.1,1,.1),
         signif(as.numeric(colnames(errs)[(1:ncol(errs))[round(seq(.1,1,.1)*ncol(errs))]]),2))
    axis(side=2,at=seq(.1,1,.1),
         signif(as.numeric(rownames(errs)[(nrow(errs):1)[round(seq(.1,1,.1)*nrow(errs))]]),2))
    text(.8,.9,
         cex=.9,
         as.expression(bquote(paste("max ",phi[T]," hit"))),
         col=rgb(1,1,1))
}
