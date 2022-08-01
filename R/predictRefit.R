predictRefit <- function(prof)
{
    require(xgboost)
    ##data("predict_refit_xgb_object_pcawg",package="ASCAT.sc")
    predict_refit_xgb_pcawg <- xgb.load(list.files(system.file('extdata',
                                                               package = 'ASCAT.sc'),
                                                   full.names = TRUE,
                                                   pattern="1.6.0.1"))
    .getTotal <- function(cn)
    {
        .getTotalDistrib <- function(cn)
        {
            nms <- as.character(0:7)
            y <- as.numeric(cn[,"total_copy_number"])
            y[y>=7] <- 7
            ny <- prop.table(tapply(1:nrow(cn),y,function(z) sum((as.numeric(cn[z,"end"])-as.numeric(cn[z,"start"]))/1000000)))
            y2 <- prop.table(tapply(1:nrow(cn),y,function(z) length(z)))
            y <- c(ny[nms],y2[nms])
            y[is.na(y)] <- 0
            names(y) <- c(paste0("sizes",0:7),paste0("numbers",0:7))
            y
        }
        cnm1 <- cn
        cnp1 <- cn
        cnm1[, "total_copy_number"] <- as.numeric(cnm1[, "total_copy_number"])-1
        cnp1[, "total_copy_number"] <- as.numeric(cnp1[, "total_copy_number"])+1
        c(.getTotalDistrib(cnm1),
          .getTotalDistrib(cn),
          .getTotalDistrib(cnp1))
    }
    predsProfs <- predict(predict_refit_xgb_pcawg,t(as.matrix(.getTotal(prof))))
}
