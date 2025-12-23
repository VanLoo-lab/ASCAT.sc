plotSolution <- function(tracksSingle,
                         purity,
                         ploidy,
                         sol=NULL,
                         ylim=c(0,8),
                         gamma=0.55,
                         ismale=F,
                         isPON=F,
                         ismedian=FALSE,
                         allchr=NULL,
                         rainbowChr=FALSE,
                         hideCN=FALSE,
                         transparentCN=FALSE,
                         zoomPoints=5,
                         ambiguousFlag=TRUE,
                         svinput=NULL,
                         is_pdf=FALSE,
                         ...)
{
  meansSeg <- fitProfile(tracksSingle,purity,ploidy,gamma=gamma, ismale=ismale, isPON=isPON, ismedian = ismedian)
  tracksSingle <- normaliseByPloidy(tracksSingle, ismedian = ismedian)
  breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))/1e+06))

  LWD_SEGS_INTEGER <- 10
  LWD_SEGS_INTEGER2 <- 5
  BASECOLOUR1 <- "grey60" ##"lightblue3"
  BASECOLOUR1 <- col2rgb(BASECOLOUR1)[,1]
  BASECOLOUR1 <- rgb(BASECOLOUR1[1]/255,BASECOLOUR1[2]/255,BASECOLOUR1[3]/255,.4)
  BASECOLOUR2 <- "black" ##"steelblue4"

  set.seed(10)
  par(plt = c(0, 1, 0, 0.93), new = TRUE, fig = c(0, 1, 0, 1))
  par(mar=c(5,5,4,1))

  plot(0, 0,
       col = rgb(0, 0, 0, 0),
       xaxt = "n",
       yaxt = "n",
       xlim = c(0, max(breaks)-80),
       xlab = "",
       ylab="",
       frame = F,
       axes=FALSE,
       ylim=ylim+c(-1,2),
       ...)

  clrs <- c("indianred1","chocolate1","orange",
            "lightgoldenrod","khaki1",  "palegreen",
            "lightgreen","seagreen1", "mediumaquamarine",
            "aquamarine","cadetblue1","turquoise2",
            "skyblue", "steelblue1","lightslateblue",
            "mediumpurple1","violet", "plum2", "plum1",
            "pink1", "lightpink","palevioletred1","lightcoral",
            "lightcoral")

  labels <- if(is.null(names(tracksSingle$lCTS))) names(breaks)[2:length(breaks)] else names(tracksSingle$lCTS)
  labels <- gsub("chr", "", labels)

  for (i in 0:ylim[2]){

    #### horizontal breaks ####

    segments(breaks[1],
             i,
             breaks[-1],
             i,
             lwd = 1,
             col = BASECOLOUR1,
             lend = 1,
             lty = 2)

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
    y_points <- ifelse(y_points>ylim[2]+.2,
                       ylim[2]+.2,
                       y_points)
    y_points <- ifelse(y_points<0,
                       0,
                       y_points)
    segments(starts/1e+06 + breaks[i] + 2,
             y_points,
             ends/1e+06 + breaks[i] + 2,
             y_points,
             col = BASECOLOUR1,
             pch = 16,
             lwd=as.numeric(zoomPoints),
             cex=1)

    nonround <- transform_bulk2tumour(sapply(meansSeg[[i]], function(x) x$mu),
                                      purity,
                                      ploidy,
                                      gamma=gamma,
                                      ismale=ismale,
                                      isPON=isPON,
                                      isX=i==23)

    segments(lSegs$loc.start[which(nonround >=0 & nonround <=ylim[2])]/1e+06 + breaks[i],
             nonround[which(nonround >=0 & nonround <=ylim[2])],
             lSegs$loc.end[which(nonround >=0 & nonround <=ylim[2])]/1e+06 + breaks[i],
             nonround[which(nonround >=0 & nonround <=ylim[2])],
             lwd = 2,
             col = rgb(0.1, 0.1, 0.1, 0.4))

    if (any(nonround < 0))
    {
      segments(lSegs$loc.start[which(nonround <0)]/1e+06 + breaks[i],
               0,
               lSegs$loc.end[which(nonround <0)]/1e+06 + breaks[i],
               0,
               lwd = 2,
               col = rgb(0.1, 0.1, 0.1, 0.4))
    }
    if(any(nonround > 8))
    {
      segments(lSegs$loc.start[which(nonround >ylim[2])]/1e+06 + breaks[i],
               ylim[2],
               lSegs$loc.end[which(nonround >ylim[2])]/1e+06 + breaks[i],
               ylim[2],
               lwd = 2,
               col = rgb(0.1, 0.1, 0.1, 0.4))
    }

    #### Main copy number segments ####

    cn <- round(sapply(meansSeg[[i]], function(x) x$roundmu))

    if(!hideCN & !transparentCN){

      if (rainbowChr){

        ##### cn between 0 and 8 #####
        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER,
                 col =clrs[i], lend=1)

        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER2,
                 col = BASECOLOUR2, lend=1)

      }
      else {

        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER,
                 col = BASECOLOUR2, lend=1)
      }

      ##### cn under 0 or over 8 #####

      if (any(cn < 0)){


        segments(lSegs$loc.start[which(cn <0)]/1e+06 + breaks[i] + 2 ,
                 0,
                 lSegs$loc.end[which(cn <0)]/1e+06 + breaks[i] + 2 ,
                 0,
                 lwd = LWD_SEGS_INTEGER,
                 col =rgb(1,0,0), lend=1)

      }
      if(any(cn > 8) ){

        segments(lSegs$loc.start[which(cn >ylim[2])]/1e+06 + breaks[i] + 2 ,
                 ylim[2],
                 lSegs$loc.end[which(cn >ylim[2])]/1e+06 + breaks[i] + 2 ,
                 ylim[2],
                 lwd = LWD_SEGS_INTEGER,
                 col =rgb(1,0,0), lend=1)

      }

    }
    else if (transparentCN & !hideCN){




      if (rainbowChr){

        ##### cn between 0 and 8 #####
        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER,
                 col = adjustcolor(clrs[i], alpha.f = 0.2), lend=1)

        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER2,
                 col =adjustcolor(BASECOLOUR2, alpha.f = 0.2), lend=1)

      }
      else {

        segments(lSegs$loc.start[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lSegs$loc.end[which(cn >=0 & cn <=ylim[2])]/1e+06 + breaks[i] + 2 ,
                 cn[which(cn >=0 & cn <=ylim[2])],
                 lwd = LWD_SEGS_INTEGER,
                 col =adjustcolor(BASECOLOUR2, alpha.f = 0.2), lend=1)
      }

      ##### cn under 0 or over 8 #####

      if (any(cn < 0)){


        segments(lSegs$loc.start[which(cn <0)]/1e+06 + breaks[i] + 2 ,
                 0,
                 lSegs$loc.end[which(cn <0)]/1e+06 + breaks[i] + 2 ,
                 0,
                 lwd = LWD_SEGS_INTEGER,
                 col =rgb(1,0,0,0.3), lend=1)

      }
      if(any(cn > 8) ){

        segments(lSegs$loc.start[which(cn >ylim[2])]/1e+06 + breaks[i] + 2 ,
                 ylim[2],
                 lSegs$loc.end[which(cn >ylim[2])]/1e+06 + breaks[i] + 2 ,
                 ylim[2],
                 lwd = LWD_SEGS_INTEGER,
                 col =rgb(1,0,0, 0.3), lend=1)

      }

    }
    #### vertical breaks ####

    if(i==1){
      segments(breaks[i],
               -1.25,
               breaks[i],
               9.25,
               lwd = 4,
               col = BASECOLOUR2, lend=1)
    }
    else{
      segments(breaks[i],
               -0.8,
               breaks[i],
               8.8,
               lwd = 4,
               col = BASECOLOUR2, lend=1)
    }

    #### chr rectangles ####

    segments(breaks[i] +2,
             9,
             breaks[i+1]+2,
             9,
             lwd = 25,
             col = BASECOLOUR2, lend=1)
    segments(breaks[i]+2,
             -1,
             breaks[i+1]+2,
             -1,
             lwd = 25,
             col = BASECOLOUR2, lend=1)

    #### chr text ####
    text(x = breaks[i] + (breaks[i+1] - breaks[i])/2,
         y = ylim[2]+1,
         labels=labels[i],
         col="white",
         cex = 1.2
         )
    text(x = breaks[i] + (breaks[i+1] - breaks[i])/2,
         y = -1,
         labels=labels[i],
         col="white",
         cex = 1.2)

    if(!is.null(svinput))
    {
        svcond <- gsub("chr","",svinput[,1])==gsub("chr","",allchr[i])
        svtmp <- svinput[svcond,]
        if(sum(svcond)>0)
            segments(breaks[i] + svtmp[,2]/1e+06,
                     0,
                     breaks[i] + svtmp[,2]/1e+06,
                     ylim[2]+.3,
                     lwd=LWD_SEGS_INTEGER/3,
                     lty=1,
                     col=rgb(1,0,0,.5))
    }
  }

  segments(tracksSingle$lSegs[[length(tracksSingle$lSegs)]]$output$loc.end[length(tracksSingle$lSegs[[length(tracksSingle$lSegs)]]$output$loc.end)]/1e+06 + breaks[i] + 4,
           -1.25,
           tracksSingle$lSegs[[length(tracksSingle$lSegs)]]$output$loc.end[length(tracksSingle$lSegs[[length(tracksSingle$lSegs)]]$output$loc.end)]/1e+06 + breaks[i] + 4,
           ylim[2]+1.25,
           lwd = 4,
           col = BASECOLOUR2, lend=1)

  ambiguous = ""

  if(!isGoodSolution(meansSeg) & is.null(sol)){ambiguous = TRUE}
  else if (isGoodSolution(meansSeg) & is.null(sol)){ambiguous = FALSE}
  else if (!is.null(sol) & !is.null(sol$ambiguous)) {ambiguous = sol$ambiguous}

  dpb <- median(unlist(lapply(tracksSingle$lCTS,function(x) x$records)),na.rm=T)
  dpb <- if(all(tracksSingle$lCTS[[1]]$records==tracksSingle$lCTS[[1]]$smoothed)) NA else dpb

  if(ambiguousFlag) {

      text(x = breaks[9],
           y=ylim[2]+ifelse(is_pdf,2.4,1.8),
           labels= paste0("purity=",
                          signif(purity,2),
                          "; average ploidy=",
                          signif(ploidy,2),
                          "; tumor ploidy=",
                          signif(getTumourPhi(ploidy,purity),2),
                          "; ambiguous=", ambiguous,
                          "; dpb=", dpb

    ), adj=0, pos=1)

  }
  else {

      text(x = breaks[9],
           y=ylim[2]+ifelse(is_pdf,2.4,1.8),
           labels= paste0("purity=",
                          signif(purity,2),
                          "; average ploidy=",
                          signif(ploidy,2),
                          "; tumor ploidy=",
                          signif(getTumourPhi(ploidy,purity),2),
                          "; large deep deletion fraction=", ambiguous,
                          "; dpb=", dpb

    ), adj=0, pos=1)

  }

  text(x = -40,
       y = c(0:ylim[2]),
       labels=c(0:ylim[2]),
       col=BASECOLOUR2,
       cex = 1.5)

  text(x = ifelse(is_pdf,-90,-70),
       y = ifelse(is_pdf,0.9,2.7),
       labels="Total copy number",
       col= BASECOLOUR2,
       cex = 1.5, adj=0, srt=90)

}
