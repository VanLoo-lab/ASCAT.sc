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
                             ord1 <- order(mat[isLocalOptima==1],decreasing=F)
                             ao1  <- ao <- ao[ord1,]
                             errs1 <- mat[ao]
                             aopp <- cbind(as.numeric(rownames(mat)[ao[,1]])*5,
                                           as.numeric(colnames(mat)[ao[,2]]))
                             clusts1 <- cutree(hclust(dist(aopp),met="ward.D2"),h=.1)
                             bests1 <- tapply(1:length(clusts1),clusts1,function(x) x[which.min(errs1[x])])
                             ao <- ao[bests1,]
                             ord <- order(mat[isLocalOptima==1][ord1][bests1],decreasing=F)
                             ao <- ao[ord,]
                             aopp <- cbind(as.numeric(rownames(mat)[ao[,1]])*5,
                                           as.numeric(colnames(mat)[ao[,2]]))
                             clusts <- cutree(hclust(dist(aopp),met="ward.D2"),h=.15)
                             errs <- mat[ao]
                             bests <- tapply(1:length(clusts),clusts,function(x) x[which.min(errs[x])])
                             bests <- bests[order(mat[ao[bests,]],decreasing=F)]
                             if(length(bests)<N) bests <- c(rep(bests[1],N-length(bests)),bests)
                             ao <- ao[bests,]
                             list(bao=ao[1:N,],
                                  ao=ao1,
                                  clusts=clusts1[ord1])  },
             error=function(e) {
                 print(e)
             },
             warning=function(w){
                 print(w)
             })
}
