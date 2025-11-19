get_fit_features <- function(fit_cn)
{
    ## author: Toby Baker, PhD, UCLA/Crick
    ## Calculate Segment Width
    fit_cn$width <- fit_cn$end - fit_cn$start
    fit_features <- list()
    for (total_copy_number in 0:6)
    {
        ## Proportion of Genome with specific total CN
        fit_features[[paste0("Prop_Genome_total_copy_number_", total_copy_number)]] <-
            weighted.mean(fit_cn$total_copy_number == total_copy_number, w = fit_cn$width,na.rm=T)
        ## Proportion of Segments with specific total CN
    }
    ## Proportion of Genome with Negative Total CN
    fit_features[["Prop_Genome_total_copy_number_Neg"]] <-
        weighted.mean(fit_cn$total_copy_number < 0, w = fit_cn$width,na.rm=T)
    ## Proportion of Genome with High Total CN
    fit_features[["Prop_Genome_total_copy_number_High"]] <-
        weighted.mean(fit_cn$total_copy_number > 6, w = fit_cn$width,na.rm=T)
    ## Calculate deviation
    deviation <- abs(fit_cn$total_copy_number_logr - fit_cn$total_copy_number)
    ## Deviation features
    fit_features[["Deviation_Med"]] <- median(deviation,na.rm=T)
    fit_features[["Deviation_Q_Low"]] <- as.numeric(quantile(deviation, 0.10,na.rm=T))
    fit_features[["Deviation_Q_High"]] <- as.numeric(quantile(deviation, 0.90,na.rm=T))
    df <- data.frame(fit_features, stringsAsFactors = FALSE)
    return(df)
}
