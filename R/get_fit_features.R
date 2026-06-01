get_fit_features <- function(fit_cn)
{
    ## author: Toby Baker, PhD, UCLA/Crick
    weighted_quantile <- function(x, w, probs, na.rm = TRUE) {
        if (na.rm) {
            keep <- !is.na(x) & !is.na(w)
            x <- x[keep]
            w <- w[keep]
        }
        if (length(x) == 0 || sum(w) == 0) {
            return(rep(NA_real_, length(probs)))
        }
        if (any(w < 0)) stop("weighted_quantile: negative weights")
        ord <- order(x)
        x <- x[ord]
        w <- w[ord]
        cw <- cumsum(w) / sum(w)
        approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
    }
    ## Calculate Segment Width
    fit_cn$width <- fit_cn$end - fit_cn$start
    fit_features <- list()
    for (total_copy_number in 0:6) {
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
    ## Width-weighted quantiles of deviation
    probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
    wq <- weighted_quantile(deviation, w = fit_cn$width,
                            probs = probs, na.rm = TRUE)
    fit_features[["Deviation_Med"]]     <- wq[3]
    fit_features[["Deviation_Q_VLow"]]  <- wq[1]
    fit_features[["Deviation_Q_Low"]]   <- wq[2]
    fit_features[["Deviation_Q_High"]]  <- wq[4]
    fit_features[["Deviation_Q_VHigh"]] <- wq[5]
    df <- data.frame(fit_features, stringsAsFactors = FALSE)
    return(df)
}
