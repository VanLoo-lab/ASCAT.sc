plotSunrise <- function(solution, localMinima = FALSE, plotClust=FALSE, is_sc=FALSE, N=10)
{
  tryCatch({
    
    plot(0, 0, col = rgb(0, 0, 0, 0), xlab = "ploidy", ylab = "purity",
         xaxt = "n", yaxt = "n", frame = F, xlim = c(0, 1), ylim = c(0,
                                                                     1))
    errs <- solution$errs
    
    errs <- errs - min(errs)
    errs.max <- max(solution$errs[!is.infinite(solution$errs)])
    errs[is.infinite(errs)] <- errs.max
    errs <- errs/errs.max
    errs <- errs[rev(seq_len(nrow(errs))), ]
    
    
    .getCol <- function(x) {
      vec <- c(sapply(seq(0, 1, length.out = 100), function(x) rgb(0.85,
                                                                   0.33, 0.33, 1 - x)), sapply(seq(0, 1, length.out = 100),
                                                                                               function(x) rgb(0.18, 0.34, 0.6, x)))
      vec[pmax(1, pmin(200, round(x * 200)))]
      
    }
    
    im <- matrix(.getCol(as.vector(1 - errs)), dim(errs)[1],
                 dim(errs)[2])
    rasterImage(as.raster(im), 0, 0, 1, 1)
    if(!is_sc) {
      
      points(which(colnames(errs) == solution$ploidy)/ncol(errs),
             1 - which(rownames(errs) == solution$purity)/nrow(errs),
             col = "chartreuse", pch = "X",cex = 1.5)
      
      
      if (localMinima){
        bao1 <- findLocalMinima(errs, N=N)
        
        ao <- bao1$bao
        
        text(ao[,2]/ncol(errs),
             1-ao[,1]/nrow(errs), label=1:nrow(ao),
             col = "white", pch = 19)
        ao <- bao1$ao
        if(plotClust){
          text(ao[,2]/ncol(errs),
               1-ao[,1]/nrow(errs), label=1:nrow(ao),
               col = RColorBrewer::brewer.pal(12,"Paired")[bao1$clusts], cex=.6)
          
        }
        
      }
    }
    
    else {
      
      points(which(colnames(errs) == solution$ploidy)/ncol(errs),
             (which(rownames(errs) == solution$purity) - 1)/(nrow(errs)-1),
             col = "chartreuse", pch = "X",cex = 1.5)
      
      
      if (localMinima){
        bao1 <- findLocalMinima(errs, N=N)
        
        ao <- bao1$bao
        
        text(ao[,2]/ncol(errs),
             (ao[,1] - 1)/(nrow(errs)-1), label=1:nrow(ao),
             col = "white", pch = 19)
        ao <- bao1$ao
        if(plotClust){
          text(ao[,2]/ncol(errs),
               1-ao[,1]/nrow(errs), label=1:nrow(ao),
               col = RColorBrewer::brewer.pal(12,"Paired")[bao1$clusts], cex=.6)
          
        }
        
      }
    }
    
    
    axis(side = 1, at = seq(0.1, 1, 0.1), signif(as.numeric(colnames(errs)[(1:ncol(errs))[round(seq(0.1,1, 0.1) * ncol(errs))]]), 2))
    
    if(!is_sc){
      
      axis(side = 2, at = seq(0.1, 1, 0.1), signif(as.numeric(rownames(errs)[(nrow(errs):1)[round(seq(0.1,1, 0.1) * nrow(errs))]]), 2))
      
    }
    else
    {
      
      axis(
        side = 2, 
        at = seq(0.1, 1, 0.9/(nrow(errs)-1)),  
     
        labels = rev(rownames(errs))  
        
      )
      
    }
    text(0.8, 0.1, cex = 0.9, as.expression(bquote(paste("max ",
                                                         phi[T], " hit"))), col = rgb(1, 1, 1))
    
    if(localMinima) return(bao1)
    
    
  },
  error=function(e) {
    print(e)
  },
  warning=function(w){
    print(w)
  }
  )
}