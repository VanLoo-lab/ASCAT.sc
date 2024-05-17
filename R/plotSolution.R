plotSolution <- function(tracksSingle,
                           purity,
                           ploidy,
                           sol=NULL,
                           ylim=c(-1,10),
                           gamma=0.55,
                           ismale=F,
                           isPON=F,
                           allchr=NULL,
                           rainbowChr=TRUE,
                           ...)
{
  meansSeg <- fitProfile(tracksSingle,purity,ploidy,gamma=gamma, ismale=ismale, isPON=isPON)
  tracksSingle <- normaliseByPloidy(tracksSingle)
  breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))
  set.seed(10)
  par(plt = c(0, 1, 0, 0.92), new = TRUE, fig = c(0, 1, 0, 1))
  
  plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
       xlim = c(0, max(breaks)-80), xlab = "",ylab="",frame = F, axes=FALSE, ylim=c(-1,10),...)
  
  
  clrs <- c("indianred1","chocolate1","orange", "lightgoldenrod","khaki1",  "palegreen", "lightgreen","seagreen1", "mediumaquamarine", "aquamarine","cadetblue1","turquoise2",
            "skyblue", "steelblue1","lightslateblue", "mediumpurple1","violet", "plum2", "plum1", "pink1", "lightpink","palevioletred1","lightcoral", "lightcoral")
  if(is.null(allchr))
    labels <- if(is.null(names(tracksSingle$lCTS))) names(breaks)[2:length(breaks)]
  else names(tracksSingle$lCTS)
  
  for (i in 0:8){
    
    #### horizontal breaks ####
    
    segments(breaks[1],
             i,
             breaks[-1],
             i,
             lwd = 1,
             col ="lightblue3", lend=1, lty = 2)
    
  }
  
  for (i in 1:length(tracksSingle$lSegs)) {
    lSegs <- tracksSingle$lSegs[[i]]$output
    lCTS <- tracksSingle$lCTS[[i]]
    
    starts <- lCTS$start
    smooth <- lCTS$smoothed
    ends <- lCTS$end
    #### Data points ####
    
    
    y_points <- transform_bulk2tumour(smooth,
                                      purity,
                                      ploidy,
                                      gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23)
    y_points <- ifelse(y_points>8.2,8.2,y_points)
    y_points <- ifelse(y_points<0,0,y_points)
    
    segments(starts/1e+06 + breaks[i] + 2,
             y_points,
             ends/1e+06 + breaks[i] + 2,
             y_points,
             
             col = rgb(0.6, 0.8, 1, 0.4),
             pch = 16, cex=1)
    
    nonround <- transform_bulk2tumour(sapply(meansSeg[[i]], function(x) x$mu),
                                      purity,
                                      ploidy,
                                      gamma=gamma, ismale=ismale, isPON=isPON, isX=i==23)
    
    if (all(nonround >= 0) & all(nonround <= 8)){
      
      segments(lSegs$loc.start/1e+06 + breaks[i] ,
               nonround,
               lSegs$loc.end/1e+06 + breaks[i],
               nonround,
               lwd = 2,
               col = rgb(0.1, 0.1, 0.1, 0.4)
               
      )
    }
    #### Main copy number segments ####
    
    cn <- round(sapply(meansSeg[[i]], function(x) x$roundmu))
    
    if (all(cn <= 8) & all(cn >= 0)){
      
      if (rainbowChr){
        
        segments(lSegs$loc.start/1e+06 + breaks[i] + 2 ,
                 cn,
                 lSegs$loc.end/1e+06 + breaks[i] + 2 ,
                 cn,
                 lwd = 17,
                 col =clrs[i], lend=1)
        
        segments(lSegs$loc.start/1e+06 + breaks[i] + 2,
                 cn,
                 lSegs$loc.end/1e+06 + breaks[i] + 2 ,
                 cn,
                 lwd = 7,
                 col ="steelblue4", lend=1)
        
      }
      
      else {
        
        segments(lSegs$loc.start/1e+06 + breaks[i] + 2,
                 cn,
                 lSegs$loc.end/1e+06 + breaks[i] + 2 ,
                 cn,
                 lwd = 17,
                 col ="steelblue4", lend=1)
      }
      
    }
    
    #### vertical breaks ####
    
    if(i==1){
      segments(breaks[i],
               -1.25,
               breaks[i],
               9.25,
               lwd = 4,
               col ="steelblue4", lend=1)
    }
    else{
      segments(breaks[i],
               -0.8,
               breaks[i],
               8.8,
               lwd = 4,
               col ="steelblue4", lend=1)
    }
    
    #### chr rectangles ####
    
    segments(breaks[i] +2,
             9,
             breaks[i+1]+2,
             9,
             lwd = 25,
             col ="steelblue4", lend=1)
    segments(breaks[i]+2,
             -1,
             breaks[i+1]+2,
             -1,
             lwd = 25,
             col ="steelblue4", lend=1)
    
    #### chr text ####
    
    if(i>1){
      
      text(x = breaks[i] + (breaks[i+1] - breaks[i])/2,
           y = 9,
           labels=labels[i],
           col="white",
           cex = 1.2
      )
      text(x = breaks[i] + (breaks[i+1] - breaks[i])/2,
           y = -1,
           labels=labels[i],
           col="white",
           cex = 1.2)
    }
    else{
      text(x = breaks[i+1]/2,
           y = 9,
           labels=labels[i],
           col="white",
           cex = 1.2
      )
      text(x = breaks[i+1]/2,
           y = -1,
           labels=labels[i],
           col="white",
           cex = 1.2)
    }
    
  }
  
  segments(tracksSingle$lSegs[[24]]$output$loc.end[2]/1e+06 + breaks[i] + 4,
           -1.25,
           tracksSingle$lSegs[[24]]$output$loc.end[2]/1e+06 + breaks[i] + 4,
           9.25,
           lwd = 4,
           col ="steelblue4", lend=1)
  
  dpb <- median(unlist(lapply(tracksSingle$lCTS,function(x) x$records)),na.rm=T)
  dpb <- if(all(tracksSingle$lCTS[[1]]$records==tracksSingle$lCTS[[1]]$smoothed)) NA else dpb
  text(x = breaks[9], y=9.8, labels= paste0("purity=",
                                            signif(purity,2),
                                            "; average ploidy=",
                                            signif(ploidy,2),
                                            "; tumor ploidy=",
                                            signif(getTumourPhi(ploidy,purity),2),
                                            if(!is.null(sol)) paste0("; ",ifelse(is.null(sol$ambiguous),"",paste0("ambiguous:",sol$ambiguous)))
                                            
  ), adj=0, pos=1)
  text(x = -40,
       y = c(0:8),
       labels=c(0:8),
       col="steelblue4",
       cex = 1.5
  )
  text(x = -90,
       y = 2.7,
       labels="Total copy number",
       col="steelblue4",
       cex = 1.5, adj=0, srt=90
  )
  
}
