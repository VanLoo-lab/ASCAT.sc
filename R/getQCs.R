getQCs <- function(res)
{
    distance_to_integer_total <- function(prof)
    {
        widths <- (prof[,"end"]-prof[,"start"])/1000000
        tot <- prof[,"total_copy_number_logr"]
        rem <- !is.na(widths) & !is.na(tot)
        dists <- sum((widths*(tot-round(tot))^2)[rem])/sum(widths[rem])
    }
    get_best_rho <- function(prof)
    {
	getBAF_dist <- function(prof, BAFs)
	{
            widths <- (prof[,"endpos"]-prof[,"startpos"])/1000000
            NAs <- round(prof[,"total_copy_number"]*BAFs)
            NBs <- round(prof[,"total_copy_number"]-NAs)
            expected <- NAs/(NAs+NBs)
            rem <- !is.na(widths) & !is.na(BAFs) & !is.na(expected)
            dists <- sqrt(sum((widths*(expected-BAFs)^2)[rem])/sum(widths[rem]))
            return(dists)
	}
	rhos <- gsub("BAF_rho","",colnames(prof)[grepl("BAF_rho",colnames(prof))])
	dists <- sapply(rhos,function(rho) getBAF_dist(prof,prof[,grep(paste0("BAF_rho",rho),colnames(prof))[1]]))
	return(as.numeric(rhos[which.min(dists)]))
    }
    distance_to_integer_AS <- function(prof)
    {
        widths <- (prof[,"endpos"]-prof[,"startpos"])/1000000
        if(any(grepl("nA_best",colnames(prof))))
        {
            NAs <- prof[,"nA_best_overdispersion"]
            NBs <- prof[,"nB_best_overdispersion"]
            best_rho <- get_best_rho(prof)
            BAFs <- prof[,paste0("BAF_rho",best_rho)]
        }
        else
        {
            NAs <- round(prof[,"total_copy_number"]*prof[,"BAF"])
            NBs <- round(prof[,"total_copy_number"]-NAs)
            BAFs <- prof[,"BAF"]
        }
        expected <- NAs/(NAs+NBs)
        rem <- !is.na(widths) & !is.na(expected)
        dists <- sqrt(sum((widths*(expected-BAFs)^2)[rem])/sum(widths[rem]))
        return(dists)
    }
    dists_AS <- sapply(res$allProfiles_AS,distance_to_integer_AS)
    res$QC_metrics <- data.frame(distance_integer_logR = sapply(res$allProfiles, distance_to_integer_total),
                                 distance_integer_BAF = if(length(dists_AS)==0) rep(NA,length(res$allProfiles)) else dists_AS)
    if("filters_data_frame"%in%names(res))
        res$QC_metrics <- cbind(res$QC_metrics, res$filters_data_frame)
    res
}
