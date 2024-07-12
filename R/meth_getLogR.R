meth_getLogR <- function(totalintensity,
                         totalintensityNormalM,
                         totalintensityNormalF,
                         sex,
                         annot,
                         GAMMA,
                         MC.CORES)
{
    getPredicted <- function(TTI, totalintensityNormal, notalreadyinpanel)
    {
        predicted <- lm(y~.-1,
                        data=data.frame(y=TTI,
                                        X=totalintensityNormal[,notalreadyinpanel]))$fitted.values
        predicted[predicted<1] <- 1
        logr <- log2(TTI/predicted)
        meds <- stats::runmed(logr,k=5001)
        keep <- abs(meds-median(meds))<0.3 & abs(logr-meds)<mad(logr-meds)*1.5
        model <- lm(y~.-1,
                        data=data.frame(y=TTI[keep],
                                        X=totalintensityNormal[keep,notalreadyinpanel]))
        predicted <- predict(model,newdata=data.frame(y=TTI,
                                                   X=totalintensityNormal[,notalreadyinpanel]))
        predicted[predicted<1] <- 1
        logr <- log2(TTI/predicted)
        logr
    }
    getPredicted <- function(TTI, totalintensityNormal, notalreadyinpanel)
    {
        predicted <- lm(y~.-1,
                        data=data.frame(y=log2(TTI),
                                        X=log2(totalintensityNormal[,notalreadyinpanel])))$fitted.values
        predicted[predicted<0] <- 0
        logr <- log2(TTI)-predicted
        meds <- stats::runmed(logr,k=5001)
        keep <- abs(meds-median(meds))<0.3 & abs(logr-meds)<mad(logr-meds)*1.5
        model <- lm(y~.-1,
                        data=data.frame(y=log2(TTI[keep]),
                                        X=log2(totalintensityNormal[keep,notalreadyinpanel])))
        predicted <- predict(model,newdata=data.frame(y=log2(TTI),
                                                   X=log2(totalintensityNormal[,notalreadyinpanel])))
        predicted[predicted<0] <- 0
        logr <- log2(TTI)-predicted
        logr
    }
    logr <- do.call("cbind",mclapply(1:ncol(totalintensity),function(x)
    {
        if(sex[x]=="male") totalintensityNormal <- totalintensityNormalM
        if(sex[x]=="female") totalintensityNormal <- totalintensityNormalF
        totalintensityNormal <- cbind(totalintensityNormal,
                                      matrixStats::rowMedians(as.matrix(totalintensityNormal)))
        notalreadyinpanel <- !colnames(totalintensityNormal)%in%colnames(totalintensity)[x]
        TTI <- totalintensity[,x]
        cat(".")
        logr <- getPredicted(TTI, totalintensityNormal,notalreadyinpanel)
        if(sex[x]=="male")
        {
            logr[annot[,1]=="X"] <- logr[annot[,1]=="X"]-GAMMA#median(logr[annot[,1]%in%c(1:22)])-GAMMA
            logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-GAMMA#median(logr[annot[,1]%in%c(1:22)])-GAMMA
        }
        if(sex[x]=="female")
            logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-3.718385
        logr
    },mc.cores=MC.CORES))
    cat("\n")
    colnames(logr) <- colnames(totalintensity); gc();
    logr
}
