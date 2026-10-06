#######################################################################
######################### RESEGMENT TOOLS #############################
#######################################################################
# Lets the app re-run ONLY the part of the pipeline that happens
# *after* the raw (PoN-normalised) logR track has been derived, using a
# different value for a chosen parameter, and preview the resulting
# profile for a single sample without redoing idat reading / PoN
# normalisation (which is by far the most expensive part).
#
# This first implementation covers `segmentation_alpha`, and is written
# for methylation-array results produced by run_methylation_array()
# (ASCAT.ma / ASCAT.scDataMeth). The functions below intentionally take
# the same names/semantics as the corresponding steps inside
# run_methylation_array() so behaviour stays identical to a full rerun
# with the new parameter value -- only cheaper, since it starts from
# the already-computed `res$logr`.
#
# Everything here is read-only with respect to the `res` object passed
# in: resegment_methyl_sample() returns a fresh track/solution/profile
# for one sample; it is up to the caller (the Shiny app) to decide
# whether/how to fold the result back into the live object.
#######################################################################

#######################################################################
## resolve_ismale
## Robustly resolves whether a given sample is male from res$sex, which
## in practice comes in several shapes depending on how the object was
## produced:
##   - one entry per sample (the common case): index normally
##   - a single value covering the WHOLE object -- common for single-
##     cell data, where every "sample" is really a cell from the same
##     patient, so sex is recorded once rather than once per cell
##   - occasionally a value that is itself a longer vector (e.g. an
##     accidentally-nested per-sample vector) -- extracted defensively
##     rather than passed straight into a comparison, since comparing a
##     length>1 vector to "male" can't be reduced to a single TRUE/FALSE
##     (and errors outright under R's stricter `&&`/`if` semantics)
## Falls back to FALSE (matching this app's existing convention
## elsewhere) whenever sex can't be resolved unambiguously, rather than
## erroring or silently mis-indexing.
#######################################################################

resolve_ismale <- function(sex, sample_index, sample_name = NULL)
{
  if (is.null(sex) || length(sex) == 0) return(FALSE)

  val <- NULL
  if (!is.null(sample_name) && !is.null(names(sex)) && sample_name %in% names(sex)) {
    val <- sex[[sample_name]]
  } else if (length(sex) == 1) {
    # Single value for the whole object -- broadcast to every sample,
    # e.g. one sex recorded for all cells of a single-cell dataset.
    val <- sex[[1]]
  } else if (!is.na(sample_index) && sample_index >= 1 && sample_index <= length(sex)) {
    val <- sex[[sample_index]]
  } else {
    return(FALSE)
  }

  # Defensive: unwrap a further nested/longer vector rather than letting
  # a length>1 comparison error out downstream.
  if (length(val) != 1) {
    if (!is.null(sample_name) && !is.null(names(val)) && sample_name %in% names(val)) {
      val <- val[[sample_name]]
    } else if (!is.na(sample_index) && sample_index >= 1 && sample_index <= length(val)) {
      val <- val[[sample_index]]
    } else if (length(val) >= 1) {
      val <- val[[1]]
    } else {
      return(FALSE)
    }
  }

  isTRUE(tolower(as.character(val)) == "male")
}

#######################################################################
## annotate_solution_with_params / append_refit_log
## Small helpers used whenever a re-segmented/refitted result is
## "Apply"-ed: attach the exact parameter values used directly onto the
## solution object (so anyone inspecting
## allSolutions.refitted.manual[[i]] later can see what produced it),
## and keep a running per-object audit trail of every such change --
## a single global field can't represent different samples having been
## refit with different parameter values.
#######################################################################

annotate_solution_with_params <- function(solution, params) {
  for (nm in names(params)) solution[[nm]] <- params[[nm]]
  solution
}

append_refit_log <- function(log, sample, tab, params) {
  entry <- c(list(sample = sample, tab = tab,
                  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             params)
  c(log, list(entry))
}

#######################################################################
## get_SBDRY_for_alpha
## Reproduces the segmentation_alpha -> SBDRY lookup used inside
## run_methylation_array(): reuse ASCAT.sc's precomputed permutation
## boundaries when available for that alpha, otherwise compute them
## on the fly via DNAcopy::getbdry().
#######################################################################

get_SBDRY_for_alpha <- function(segmentation_alpha, nperms = 10000)
{
  if (!exists("SBDRYs", inherits = TRUE))
    data("SBDRYs_precomputed", package = "ASCAT.sc")

  if (as.character(segmentation_alpha) %in% names(SBDRYs))
    return(SBDRYs[[as.character(segmentation_alpha)]])

  max.ones <- floor(nperms * segmentation_alpha) + 1
  DNAcopy::getbdry(eta = 0.05, nperm = nperms, max.ones = max.ones)
}

#######################################################################
## resegment_methyl_sample
##
## res           : the ascat.sc/ascat.ma result object as held by the
##                 app (`reslocal`). Must contain, at minimum:
##                   - res$logr                : raw PoN-normalised logR
##                                                matrix (probes x samples)
##                   - res$annotations.probes  : probe annotation data
##                                                frame (chr/pos columns),
##                                                row-aligned with res$logr
##                   - res$chr                 : chromosomes used (allchr)
##                   - res$sex                 : per-sample sex ("male"/"female")
##                   - res$gamma, res$min.width: fit constants from the
##                                                original run
## sample_index  : integer index of the sample (matches
##                 names(res$allTracks.processed)/colnames(res$logr))
## segmentation_alpha : the new segmentation penalty to try
## min.width     : DNAcopy min.width; defaults to res$min.width
## purs/ploidies/maxtumourpsi : purity/ploidy grid-search ranges used by
##                 searchGrid(); these are NOT stored on the result
##                 object by run_methylation_array(), so they default to
##                 the same ranges run_methylation_array() itself
##                 defaults to. Pass explicit values if the original run
##                 used something else.
## gamma         : overrides the gamma used for fitting (searchGrid /
##                 fitProfile). Defaults to res$gamma (falling back to
##                 0.55, run_methylation_array()'s own default) when NULL.
##
## Returns list(track = <ASCAT.sc track>, solution = <searchGrid result>,
##              profile = <getProfile result>, params = list(...)) for
## that one sample. `params` records the exact parameter values this
## particular fit used (segmentation_alpha, gamma), so the caller can
## save a record of them alongside the profile.
#######################################################################

resegment_methyl_sample <- function(res,
                                     sample_index,
                                     segmentation_alpha,
                                     gamma        = NULL,
                                     min.width    = NULL,
                                     purs         = seq(0.1, 1, 0.01),
                                     ploidies     = seq(1.7, 4, 0.01),
                                     maxtumourpsi = 5,
                                     MC.CORES     = 1)
{
  if (is.null(res$logr))
    stop("res$logr (raw logR track) not found -- cannot re-segment without it")
  if (!is.matrix(res$logr) && !is.data.frame(res$logr))
    stop("res$logr must be a matrix (probes x samples, samples as columns)")

  ## Resolve the sample by NAME, not by column position: res$logr is a
  ## probes x samples matrix, and its column order is not guaranteed to
  ## match the order of res$allTracks.processed / res$sex, so we always
  ## go through the sample name shown in the app (never a bare index into
  ## the matrix itself).
  samp <- names(res$allTracks.processed)[sample_index]
  if (is.null(samp) || is.na(samp))
    stop("sample_index is out of range for res$allTracks.processed")
  if (!samp %in% colnames(res$logr))
    stop(paste0("Sample '", samp, "' is not a column of res$logr -- ",
                "check that colnames(res$logr) match the sample names ",
                "used elsewhere in the result object"))

  annot <- res$annotations.probes
  if (is.null(annot))
    stop("res$annotations.probes not found -- cannot recover probe chr/pos")

  GAMMA <- if (!is.null(gamma)) gamma else if (!is.null(res$gamma)) res$gamma else 0.55
  if (is.null(min.width))
    min.width <- if (!is.null(res$min.width)) res$min.width else 5
  allchr <- if (!is.null(res$chr)) res$chr else c(1:22, "X", "Y")

  starts <- as.numeric(as.character(annot[, "pos"]))
  ends   <- starts
  chrs   <- as.character(annot[, "chr"])

  raw_logr <- res$logr[, samp]

  ## `res$bins` is only ever populated (with an S4 object carrying a
  ## `bins` slot) when the original run used the conumee/binned branch;
  ## it stays NULL for the standard non-binned branch. Mirror that here
  ## so a binned original run gets re-segmented the same way.
  if (!is.null(res$bins))
  {
    require(ASCAT.ma)
    input <- meth_bin(raw_logr, starts = starts, ends = ends,
                       chrs = paste0("chr", chrs), res$bins@bins)
    SBDRY <- NULL  # getTrackForAll.bins recomputes/handles SBDRY internally in this branch
    track <- getTrackForAll.bins(input[[1]], input[[4]], input[[2]], input[[3]],
                                  segmentation_alpha = segmentation_alpha,
                                  min.width = min.width,
                                  allchr = gsub("chr", "", allchr))
  }
  else
  {
    require(ASCAT.scDataMeth)
    SBDRY <- get_SBDRY_for_alpha(segmentation_alpha)
    .logr <- meth_winsorise_ascat(raw_logr)
    track <- getTrackForAll.bins(.logr, paste0("chr", chrs), starts, ends,
                                  segmentation_alpha = segmentation_alpha,
                                  transform = FALSE, ismedian = TRUE,
                                  min.width = min.width, SBDRY = SBDRY,
                                  allchr = gsub("chr", "", allchr))
  }

  ismale <- resolve_ismale(res$sex, sample_index, samp)

  sol <- try(searchGrid(track, purs = purs, ploidies = ploidies,
                         maxTumourPhi = maxtumourpsi, gamma = GAMMA,
                         ismale = ismale, ismedian = TRUE), silent = TRUE)
  if (inherits(sol, "try-error"))
    stop(paste0("Could not fit purity/ploidy with segmentation_alpha=", segmentation_alpha,
                ": ", attr(sol, "condition")$message))

  profile <- try(getProfile(fitProfile(track, purity = sol$purity, ploidy = sol$ploidy,
                                        gamma = GAMMA, ismedian = TRUE, ismale = ismale),
                             CHRS = allchr), silent = TRUE)
  if (inherits(profile, "try-error"))
    stop(paste0("Could not build the profile with segmentation_alpha=", segmentation_alpha,
                ": ", attr(profile, "condition")$message))

  list(track = track, solution = sol, profile = profile,
       params = list(segmentation_alpha = segmentation_alpha, gamma = GAMMA))
}
