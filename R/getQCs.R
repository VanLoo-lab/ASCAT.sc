getQCs <- function(res)
{
    distance_to_integer_total <- function(prof)
    {
        widths <- (prof[,"end"]-prof[,"start"])/1000000
        tot <- prof[,"total_copy_number_logr"]
        rem <- !is.na(widths) & !is.na(tot)
        dists <- sum((widths*(tot-round(tot))^2)[rem])/sum(widths[rem])
    }
    mad_scaled <- function(x)
    {
      med <- median(x, na.rm = TRUE)
      1.4826 * median(abs(x - med), na.rm = TRUE) # scaled with sd
    }
    mapd <- function(x, chr)
    {
      diffs <- unlist(
        tapply(x, chr, function(v) {
          if(length(v) < 2) return(NULL)
             abs(diff(v))
        }),
        use.names = FALSE
      )
      median(diffs, na.rm = TRUE)
    }
    spikiness <- function(x, chr)
    {
      vals <- unlist(
        tapply(x, chr, function(v) {
          if(length(v) < 3) return(NULL)
          abs(diff(diff(v)))
        }),
        use.names = FALSE
      )
      median(vals, na.rm = TRUE)
    }
    autocorr_within_chr <- function(x, chr)
    {
      split_vals <- split(x, chr)
      per_chr <- sapply(split_vals, function(v) {
        if(length(v) < 3) return(NA_real_)
        cor(v[-length(v)], v[-1], use = "complete.obs")
      })
      weighted.mean(per_chr,
                    weights = sapply(split_vals, length),
                    na.rm = TRUE)
    }
    extra_qc <- lapply(seq_along(res$allProfiles), function(i)
      {
      prof <- res$allProfiles[[i]]
      lcts <- res$allTracks.processed[[i]]$lCTS
      
      bins <- data.frame()
      
      for(chr in names(lcts))
      {
        bin_df <- lcts[[chr]]
        
        seg <- prof[prof$chromosome == chr, ]
        
        bins_chr <- data.frame(
          chr = chr,
          logR = bin_df$smoothed,
          records = bin_df$records,
          state = rep(seg$total_copy_number, seg$num.mark),
          seg_logr = rep(seg$logr, seg$num.mark)
        )
        
        bins <- rbind(bins, bins_chr)
      }
      
      bins$resid_seg <- bins$logR - bins$seg_logr
      
      breaks <- c(
        FALSE, 
        diff(bins$state) != 0 &
        bins$chr[-1] == bins$chr[-nrow(bins)]
      )
      
    rho <- res$allSolutions[[i]]$purity
    psi <- res$allSolutions[[i]]$ploidy
    expected_logR <- function(n, rho, psi)
    {
      log2(
        (rho * n + 2 * (1-rho)) /
        (rho * psi + 2 * (1-rho))
      )
    }
    bins$expected <-
      expected_logR(
        bins$state,
        rho,
        psi
      )

    seg_resid_int <-
      tapply(
        bins$seg_logr - bins$expected,
        interaction(bins$chr, bins$state),
        median,
        na.rm = TRUE
      )

      data.frame(
        MAPD_gc_corrected =
          mapd(bins$logR, bins$chr),

        spikiness =
          spikiness(bins$logR, bins$chr),

        MBRSM_dispersion =
          mad_scaled(bins$resid_seg),

        breakpoints =
          sum(breaks, na.rm = TRUE),

        state_mode =
          as.integer(names(which.max(table(bins$state)))),

        autocorr =
          autocorr_within_chr(
            bins$resid_seg,
            bins$chr
          ),

        MSRSI_non_integerness =
          mad_scaled(seg_resid_int)
      )
    })
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
    dists_AS <- sapply(res$allProfiles_AS, function(x) distance_to_integer_AS(x$nprof.fixed))
    res$QC_metrics <- data.frame(distance_integer_logR = sapply(res$allProfiles, distance_to_integer_total),
                                 distance_integer_BAF = if(length(dists_AS)==0) rep(NA,length(res$allProfiles)) else dists_AS)
    extra_qc <- do.call(rbind, extra_qc)
    res$QC_metrics <- cbind(res$QC_metrics, extra_qc)
    if("filters_data_frame"%in%names(res))
        res$QC_metrics <- cbind(res$QC_metrics, res$filters_data_frame)
    res
}
