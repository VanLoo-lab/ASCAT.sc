meth_derive_pon <- function(idat_dir,
                            allchr=c(1:22,"X","Y"),
                            GAMMA=0.55)
{
    inferSex <- function(tt, annot)
    {
        chrX <- apply(tt[annot[,1]=="X",],2,function(x) median(log2(x)))
        chrY <- apply(tt[annot[,1]=="Y",],2,function(x) median(log2(x)))
        autos <- apply(tt[!annot[,1]%in%c("X","Y"),],2,function(x) median(log2(x)))
        sex <- ifelse(chrX-autos>0,ifelse(chrY-autos< -1, "female",NA), ifelse(chrY-autos< -1, NA,"male"))
        sex[is.na(sex)] <- "female" ## default sex
        sex[chrY-autos> -1] <- "male"
        sex
    }

    derivePON <- function(tt, annot, madtimes=1.5)
    {
        NN <- nrow(tt)
        chrX <- annot[,"chr"]=="X"
        chrY <- annot[,"chr"]=="Y"
        chrA <- !annot[,"chr"]%in%c("Y","X")
        NNX <- sum(chrX)
        NNA <- sum(chrA)
        NNY <- sum(chrY)
        ranksA <- rowSums(apply(tt[chrA,],2,function(x)
        {
            cat(".")
            rr <- rank(x)
            log(rr/(NNA-rr+1))
        }))
        ranksX <- rowSums(apply(tt[chrX,],2,function(x)
        {
            cat(".")
            rr <- rank(x)
            log(rr/(NNX-rr+1))
        }))
        ranksY <- rowSums(apply(tt[chrY,],2,function(x)
        {
            cat(".")
            rr <- rank(x)
            log(rr/(NNY-rr+1))
        }))
        badlociA <- abs(ranksA)>median(ranksA)+mad(ranksA)*madtimes
        badlociX <- abs(ranksX)>median(ranksX)+mad(ranksX)*madtimes
        badlociY <- abs(ranksY)>median(ranksY)+mad(ranksY)*madtimes
        badloci <- rep(F,nrow(tt))
        badloci[(which(chrA))[badlociA]] <- T
        badloci[(which(chrX))[badlociX]] <- T
        badloci[(which(chrY))[badlociY]] <- T
        list(badloci=badloci,
             ranks=c(ranksA,ranksX,ranksY))
    }


    winsorise_ascat <- function(x)
    {
        madWins <- function(x,tau,k)
        {
            medianFilter <- function(x,k){
                n <- length(x)
                filtWidth <- 2*k + 1
                                        #Make sure filtWidth does not exceed n
                if(filtWidth > n){
                    if(n==0){
                        filtWidth <- 1
                    }else if(n%%2 == 0){
                                        #runmed requires filtWidth to be odd, ensure this:
                        filtWidth <- n - 1
                    }else{
                        filtWidth <- n
                    }
                }
                runMedian <- stats::runmed(x,k=filtWidth,endrule="median")
                return(runMedian)
            }
            xhat <- medianFilter(x,k)
            d <- x-xhat
            SD <- mad(d)
            z <- tau*SD
            xwin <- xhat + psi(d, z)
            outliers <- rep(0, length(x))
            outliers[x > xwin] <- 1
            outliers[x < xwin] <- -1
            return(list(ywin=xwin,sdev=SD,outliers=outliers))
        }
        psi <- function(x,z)
        {
            xwin <- x
            xwin[x < -z] <- -z
            xwin[x > z] <- z
            return(xwin)
        }
        lrwins = vector(mode="numeric",length=length(x))
        lrwins[is.na(x)] = NA
        lrwins[!is.na(x)] = madWins(x[!is.na(x)],3,15)$ywin
        lrwins
    }
    ## ##################################################

    ## ##################################################
    suppressPackageStartupMessages(require(minfi))
    suppressPackageStartupMessages(require(conumee))
    ## ##################################################
    print("## read idat directory")
    rgSet <- read.metharray.exp(idat_dir,force=T); gc();
    ## ##################################################
    print("## preprocess raw data - extract unmeth/meth signal unprocessed")
    data <- preprocessRaw(rgSet);
    ## ##################################################
    totalintensity <- getMeth(data)+getUnmeth(data)+1; gc();
    annot <- as.data.frame(getAnnotation(data));
    annot[,"chr"] <- gsub("chr","",as.character(annot[,"chr"]))
    annot <- annot[as.character(annot[,"chr"])%in%as.character(allchr),]
    ## ##################################################
    ## order chromosome, starts and ends of probes
    annot <- annot[order(annot[,"chr"],annot[,"pos"],decreasing=F),]
    starts <- as.numeric(as.character(annot[,"pos"]))
    ends <- as.numeric(as.character(annot[,"pos"]))
    chrs <- as.character(annot[,"chr"])
    ## ##################################################
    ## Final ordered total intensities for all and normal
    totalintensity <- totalintensity[rownames(annot),]
    sex <- inferSex(totalintensity,annot)
    id_normals <- colnames(totalintensity)[sex=="female"]
    totalintensityNormalF <- totalintensity[,colnames(totalintensity)%in%id_normals]; gc();
    id_normals <- colnames(totalintensity)[sex=="male"]
    totalintensityNormalM <- totalintensity[,colnames(totalintensity)%in%id_normals]; gc();
    ## ##################################################
    ## ##################################################
    print("## derive PoN-normalised logr")
    ## log PoN-fitted intensity values of all probes for all samples
    logr <- sapply(1:ncol(totalintensity),function(x)
    {
        if(sex[x]=="male") totalintensityNormal <- totalintensityNormalM
        if(sex[x]=="female") totalintensityNormal <- totalintensityNormalF
        notalreadyinpanel <- !colnames(totalintensityNormal)%in%colnames(totalintensity)[x]
        TTI <- totalintensity[,x]
        cat(".")
        predicted <- lm(y~.-1,
                        data=data.frame(y=TTI,
                                        X=totalintensityNormal[,notalreadyinpanel]))$fitted.values
        cat("-")
        predicted[predicted<1] <- 1
        logr <- log2(totalintensity[,x]/predicted)
        if(sex[x]=="male")
        {
            logr[annot[,1]=="X"] <- logr[annot[,1]=="X"]-median(logr[annot[,1]%in%c(1:22)])-GAMMA
            logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-median(logr[annot[,1]%in%c(1:22)])-GAMMA
        }
        if(sex[x]=="female")
            logr[annot[,1]=="Y"] <- logr[annot[,1]=="Y"]-3.718385
        logr
    })
    rownames(logr) <- rownames(annot)
    colnames(logr) <- colnames(totalintensity); gc();
    ## ##################################################
    ## ##################################################
    badloci <- derivePON(logr,annot=annot)
    ## ##################################################
    list(totalintensity=totalintensity,
         badloci.annot=annot[badloci$badloci,])
}

