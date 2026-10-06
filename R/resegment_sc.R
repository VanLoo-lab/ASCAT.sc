#######################################################################
##################### SC / SHALLOW-COVERAGE REFIT TOOLS ###############
#######################################################################
# Companion to resegment_methyl.R, for results produced by
# run_sc_sequencing() (used for both true single-cell and bulk
# shallow-coverage WGS data -- both share the same object structure and
# only differ in how many "samples"/cells they contain).
#
# Unlike the methylation-array case, changing the parameters here means
# re-running TWO steps, not one: re-binning the raw per-sample coverage
# track at a new `binsize`, and re-segmenting the result with a new
# `segmentation_alpha`. Both are cheap compared to re-reading BAMs, so
# this starts from the raw (native-resolution) coverage counts already
# stored on the result object (`res$allTracks[[sample]]$lCTS.tumour`)
# rather than re-reading anything from disk.
#
# IMPORTANT CAVEAT: run_sc_sequencing() builds a panel-of-normals
# correction (`lNormals`) internally when normal BAMs are supplied, but
# does NOT keep `lNormals`/`lCTS.normal` on the final returned object --
# only `res$isPON` (a flag) survives. That means if the original run
# used a panel of normals, this preview cannot exactly reproduce that
# correction and re-bins/re-segments WITHOUT it. resegment_sc_sample()
# surfaces this via the returned `note` field so the app can warn the
# user rather than silently producing a subtly different result.
#######################################################################

#######################################################################
## .reload_reference_tracks
## Reloads the *raw*, native-resolution reference selection (`lSe`) and
## GC-content (`lGCT`) tracks for a genome build, exactly as
## run_sc_sequencing() does at the start of a run. These are NOT the
## same as res$lSe/res$lGCT stored on the result object -- those are
## already binned at the ORIGINAL binsize -- so they have to be reloaded
## from the package data to be re-binned at a new binsize.
#######################################################################

.reload_reference_tracks <- function(build, allchr, chrstring_bam = "")
{
  if (build == "hg19")
  {
    data("lSe_filtered_30000.hg19", package = "ASCAT.sc")
    data("lGCT_filtered_30000.hg19", package = "ASCAT.sc")
    allchr. <- gsub("chr", "", allchr)
    lSe  <- lapply(allchr., function(chr) lSe.hg19.filtered[[chr]])
    names(lSe) <- allchr
    names(lGCT.hg19.filtered) <- names(lSe)
    lGCT <- lapply(allchr, function(chr) lGCT.hg19.filtered[[chr]])
    START_WINDOW <- 30000
  }
  else if (build == "hg38")
  {
    data("lSe_filtered_30000.hg38", package = "ASCAT.sc")
    data("lGCT_filtered_30000.hg38", package = "ASCAT.sc")
    names(lGCT.hg38.filtered) <- names(lSe.hg38.filtered)
    allchr. <- paste0("chr", gsub(chrstring_bam, "", allchr))
    lSe  <- lapply(allchr., function(chr) lSe.hg38.filtered[[chr]])
    names(lSe) <- allchr
    lGCT <- lapply(allchr., function(chr) lGCT.hg38.filtered[[chr]])
    names(lGCT) <- allchr
    if (chrstring_bam == "")
      names(lGCT) <- names(lSe) <- gsub("chr", "", names(lSe))
    START_WINDOW <- 30000
  }
  else if (build == "mm39")
  {
    data("lSe_unfiltered_5000.mm39", package = "ASCAT.sc")
    data("lGCT_unfiltered_5000.mm39", package = "ASCAT.sc")
    names(lGCT)[1:length(allchr)] <- names(lSe)[1:length(allchr)] <- allchr
    lSe  <- lapply(allchr, function(x) lSe[[x]])
    lGCT <- lapply(allchr, function(x) lGCT[[x]])
    names(lGCT) <- names(lSe) <- allchr
    START_WINDOW <- 5000
  }
  else
  {
    stop(paste0("Unsupported genome build for re-binning: '", build,
                "' (only hg19, hg38, mm39 are currently supported)"))
  }

  list(lSe = lSe, lGCT = lGCT, START_WINDOW = START_WINDOW)
}

#######################################################################
## resegment_sc_sample
##
## res           : the ascat.sc result object as held by the app
##                 (`reslocal`), as returned by run_sc_sequencing().
##                 Must contain, at minimum:
##                   - res$allTracks[[sample]]$lCTS.tumour : raw,
##                     native-resolution per-chromosome coverage counts
##                   - res$build   : "hg19" | "hg38" | "mm39"
##                   - res$chr     : chromosomes used (allchr)
##                   - res$sex     : per-sample sex ("male"/"female")
## sample_index  : integer index of the sample (matches
##                 names(res$allTracks.processed)/names(res$allTracks))
## segmentation_alpha : the new segmentation penalty to try
## binsize       : the new bin size (bp) to re-bin the raw counts at
## purs/ploidies/maxtumourpsi : purity/ploidy grid-search ranges used by
##                 searchGrid(); default to the same ranges
##                 run_sc_sequencing() itself defaults to.
##
## Returns list(track =, solution =, profile =, note =) for that one
## sample; `note` is a caveat string (see header) or NULL. Does not
## mutate `res`.
#######################################################################

resegment_sc_sample <- function(res,
                                 sample_index,
                                 segmentation_alpha,
                                 binsize,
                                 purs         = seq(0.1, 1, 0.01),
                                 ploidies     = seq(1.7, 5, 0.01),
                                 maxtumourpsi = 5,
                                 MC.CORES     = 1)
{
  if (is.null(res$allTracks))
    stop("res$allTracks (raw per-sample coverage tracks) not found -- cannot re-bin without it")

  samp <- names(res$allTracks.processed)[sample_index]
  if (is.null(samp) || is.na(samp))
    stop("sample_index is out of range for res$allTracks.processed")
  if (!samp %in% names(res$allTracks))
    stop(paste0("Sample '", samp, "' is not present in res$allTracks"))

  raw_track <- res$allTracks[[samp]]
  if (is.null(raw_track) || is.null(raw_track$lCTS.tumour))
    stop(paste0("res$allTracks[['", samp, "']]$lCTS.tumour (raw, unbinned coverage) ",
                "not found -- cannot re-bin without it"))

  if (is.null(res$build))
    stop("res$build (genome build) not found -- needed to reload reference bin tracks")
  allchr <- if (!is.null(res$chr)) res$chr else stop("res$chr not found")
  chrstring_bam <- if (!is.null(res$chrstring_bam)) res$chrstring_bam else ""

  refs <- .reload_reference_tracks(res$build, allchr, chrstring_bam)
  START_WINDOW <- refs$START_WINDOW

  if (binsize < START_WINDOW)
    stop(paste0("binsize must be >= ", START_WINDOW, " (the native bin resolution for build '",
                res$build, "') -- got ", binsize))

  nlGCT <- treatGCT(refs$lGCT, window = ceiling(binsize / START_WINDOW))
  nlSe  <- treatlSe(refs$lSe,  window = ceiling(binsize / START_WINDOW))

  nlCTS.tumour <- treatTrack(lCTS = raw_track$lCTS.tumour,
                              window = ceiling(binsize / START_WINDOW))

  ## Panel-of-normals correction can't be exactly reproduced -- see the
  ## caveat in the file header -- so this re-bins/re-segments without it.
  note <- if (isTRUE(res$isPON))
    paste0("The original run used a panel of normals for GC/mappability correction. ",
           "That panel isn't stored on the loaded object, so this preview re-bins and ",
           "re-segments WITHOUT it, and may differ slightly from a true full rerun.")
  else NULL

  SBDRY <- get_SBDRY_for_alpha(segmentation_alpha)

  svinput_arg <- if (!is.null(res$lSVinput)) res$lSVinput[[samp]] else NULL

  track <- getTrackForAll(bamfile = NULL, window = NULL,
                           lCT = nlCTS.tumour, lSe = nlSe, lGCT = nlGCT,
                           lNormals = NULL, allchr = allchr,
                           sdNormalise = 0, SBDRY = SBDRY,
                           svinput = svinput_arg,
                           segmentation_alpha = segmentation_alpha)

  ismale <- resolve_ismale(res$sex, sample_index, samp)  # defined in resegment_methyl.R

  sol <- try(searchGrid(track, purs = purs, ploidies = ploidies,
                         maxTumourPhi = maxtumourpsi,
                         ismale = ismale, isPON = FALSE), silent = TRUE)
  if (inherits(sol, "try-error"))
    stop(paste0("Could not fit purity/ploidy at binsize=", binsize,
                ", segmentation_alpha=", segmentation_alpha, ": ",
                attr(sol, "condition")$message))

  profile <- try(getProfile(fitProfile(track, purity = sol$purity, ploidy = sol$ploidy,
                                        ismale = ismale), CHRS = allchr), silent = TRUE)
  if (inherits(profile, "try-error"))
    stop(paste0("Could not build the profile at binsize=", binsize,
                ", segmentation_alpha=", segmentation_alpha, ": ",
                attr(profile, "condition")$message))

  list(track = track, solution = sol, profile = profile, note = note,
       params = list(segmentation_alpha = segmentation_alpha, binsize = binsize))
}
