findLocalMinima <- function(mat, N=5)
{
  tryCatch({shifts <- list(c(1,1),
                           c(0,1),
                           c(1,0),
                           c(-1,-1),
                           c(-1,+1),
                           c(-1,0),
                           c(0,-1),
                           c(1,-1))
  nr <- nrow(mat)
  nc <- ncol(mat)
  applyShifts <- function(mat,shifts,nr,nc)
  {
    if(shifts[1]==-1)
      mat <- rbind(mat[-c(1),],Inf)
    if(shifts[1]==1)
      mat <- rbind(Inf,mat[-c(nr),])
    if(shifts[2]==-1)
      mat <- cbind(mat[,-c(1)],Inf)
    if(shifts[2]==1)
      mat <- cbind(Inf,mat[,-c(nc)])
    mat
  }
  nmat <- list()
  for(i in 1:8)
  {
    nmat[[i]] <- mat
    nmat[[i]] <- apply(applyShifts(nmat[[i]],shifts[[i]],nr,nc)>mat,2,as.numeric)
  }
  
  isLocalOptima <- Reduce("*",nmat)
  ao <- which(isLocalOptima==1,arr.ind=T)
  aopp <- cbind(as.numeric(rownames(mat)[ao[,1]])*2,
                as.numeric(colnames(mat)[ao[,2]]))
  clusts1 <- cutree(hclust(dist(aopp),met="ward.D2"),h=.25)
  ord1 <- order(mat[isLocalOptima==1],decreasing=F)
  ao1  <- ao <- ao[ord1,]
  ao <- ao[!duplicated(clusts1[ord1]),]
  aopp <- cbind(as.numeric(rownames(mat)[ao[,1]])*2,
                as.numeric(colnames(mat)[ao[,2]]))
  clusts <- cutree(hclust(dist(aopp),met="ward.D2"),h=.25)
  ord <- order(mat[isLocalOptima==1][ord1][!duplicated(clusts1[ord1])],decreasing=F)
  ao <- ao[ord,]
  list(bao=ao[!duplicated(clusts[ord]),][1:N,],
       ao=ao1,
       clusts=clusts1[ord1])  },
  error=function(e) {
    print(e)
  },
  warning=function(w){
    print(w)
  })
  
}