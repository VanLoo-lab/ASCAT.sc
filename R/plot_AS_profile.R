plot_AS_profile  <-  function(t,
                              col1=rgb(47/255,79/255,79/255,1),##col of minor allele
                              col2=rgb(230/255,159/255,0,1),##col of total alleles
                              lwdSegs=5.5) ##width of the segments
{
    t$chr <- as.character(t$chr)
    add <- cumsum(sapply(unique(t$chr),function(x)
    {
        max(t$endpos[t$chr==x])/1000
    }))
    add <- c(0,add)
    names(add) <- c(unique(t$chr),"NO")
    notna <- function(x) ifelse(is.na(x),0,x)
    t$averageMin <- t$nB
    t$averageMaj <- t$nA
    t$averageTot <- t$averageMin+t$averageMaj
    t$averageTot[t$averageTot>10] <- 10
    par(mar=c(1,3,1,1))
    plot(0,0,xaxt="n",yaxt="n",
         xlab=paste0(""),
         ylab="",
         main=paste0(""),
         col=rgb(0,0,0,0),
         frame=T,
         xlim=c(0,add[length(add)]),
         ylim=c(0,10))
    title(ylab="Copy Number",line=2)
    for(i in unique(t$chr))
    {
        subt <- t[t$chr==i,]
        segments(subt$startpos/1000+add[i],
                 subt$averageMin,
                 subt$endpos/1000+add[i],
                 subt$averageMin,col=col1,lwd=lwdSegs, lend="butt")
        segments(subt$startpos/1000+add[i],
                 subt$fitted,
                 subt$endpos/1000+add[i],
                 subt$fitted,col=rgb(.1,.1,.1,.1),lwd=lwdSegs*.8, lend="butt")
        segments(subt$startpos/1000+add[i],
                 subt$averageTot,
                 subt$endpos/1000+add[i],
                 subt$averageTot,col=col2,lwd=lwdSegs, lend="butt")
    }
    abline(v=add,lty=1,col=rgb(0,0,0,.5))
    text(unique(t$chr),
         x=add[-c(1)]-25000,
         y=9.1,
         cex=.6)
    axis(side=2,at=0:12, las=2)
    abline(h=0:10,lty=1,col=rgb(.6,.6,.6,.3))
    return(add)
}
