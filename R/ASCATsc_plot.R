# =============================================================================
# ASCAT.sc Unified Plotting Function
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 1: Load Required Packages
# -----------------------------------------------------------------------------

# Core packages
required_packages <- c("data.table", "ggplot2", 
                       "cowplot", "dplyr", "ggdendro", "patchwork", "scales", "RColorBrewer", "ggsci")

# Check and load packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.\n
                  Install it with: install.packages('%s')", pkg, pkg))
  }
}

# Load packages explicitly
library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggdendro)
library(patchwork)
library(scales)
library(ggsci)

# Set default theme
theme_set(theme_cowplot())

# -----------------------------------------------------------------------------
# STEP 2: Define Plotting Function
# -----------------------------------------------------------------------------

# Get non-allele-specific colors (total CN only)
get_total_cn_colors <- function() {
  colors <- c(
    "#153570",  # 0 
    "#577aba",  # 1 
    "#c1c1c1",  # 2 
    "#e3b55f",  # 3 
    "#d6804f",  # 4 
    "#b3402e",  # 5 
    "#821010",  # 6 
    "#6a0936",  # 7 
    "#ab1964",  # 8 
    "#b6519f",  # 9 
    "#ad80b9",  # 10 
    "#c2a9d1"   # 10+ 
  )
  names(colors) <- c(as.character(0:10), "10+")
  return(colors)
}

# Get allele-specific colors
# Gradient within each total CN group
get_allele_specific_colors <- function() {
  colors <- c(
    # Total CN = 0
    "0|0" = "#153570",
    
    # Total CN = 1 (teal/cyan family - light to dark)
    "1|0" = "#345d8e",
    
    # Total CN = 2 (gray - light to dark)
    "1|1" = "#c1c1c1",
    "2|0" = "#737373",
    
    # Total CN = 3 (gold/amber - light to dark)
    "2|1" = "#edd78e",
    "3|0" = "#c49631",
    
    # Total CN = 4 (green - light to dark)
    "2|2" = "#a8d9a8",
    "3|1" = "#5fa85f",
    "4|0" = "#2d7a2d",
    
    # Total CN = 5 (coral/salmon - light to dark)
    "3|2" = "#f2b4a8",
    "4|1" = "#cf3a3a",
    "5|0" = "#7d1111",
    
    # Total CN = 6 (brown/tan - light to dark)
    "3|3" = "#d4b896",
    "4|2" = "#a67c52",
    "5|1" = "#7a5230",
    "6|0" = "#4a2c12",
    
    # Total CN = 7 (purple - light to dark)
    "4|3" = "#c7aed5",
    "5|2" = "#9b6bb5",
    "6|1" = "#6e3a8e",
    "7|0" = "#4a1570",
    
    # Overflow
    "7+" = "#3d0d5c",
    "NA" = "#BEBEBE"
  )
  return(colors)
}

# -----------------------------------------------------------------------------
# Step 3: Helper functions for states
# -----------------------------------------------------------------------------

#' Encode major|minor alleles to a display state
#' @param major Vector of major allele copy numbers
#' @param minor Vector of minor allele copy numbers
#' @return Character vector of states like "2|1", "3|0", "NA", etc.
encode_allele_state <- function(major, minor) {
  # Handle NA values first
  is_na <- is.na(major) | is.na(minor)
  
  # Ensure major >= minor (by definition, major is the larger)
  maj <- pmax(major, minor, na.rm = FALSE)
  min <- pmin(major, minor, na.rm = FALSE)
  
  total <- maj + min
  
  # Build state string directly
  state <- paste0(maj, "|", min)
  
  # Anything with total > 7 becomes "7+"
  state[!is_na & total > 7] <- "7+"
  
  # Handle NAs
  state[is_na] <- "NA"
  
  return(state)
}

# -----------------------------------------------------------------------------
# Step 3b: Map RNA cell types to DNA barcodes & grouped clustering
# -----------------------------------------------------------------------------

#' Map RNA (GEX) cell type annotations to DNA (ATAC) barcodes
#' Uses a barcode mapping file (ATAC,GEX csv) from 10X Genomics whitelist and a cell type annotation file
#'
#' @param cell_type_file Path to tab-delimited file with columns "barcode" (RNA) and "cell_type"
#' @param barcode_map_file Path to CSV with columns "ATAC" (DNA) and "GEX" (RNA)
#' @return data.table with columns "barcode" (DNA barcodes) and "cell_type"
map_cell_types_to_dna <- function(cell_type_file, barcode_map_file) {
  # Read cell type annotations (RNA barcodes)
  ct <- fread(cell_type_file, header = TRUE)
  setnames(ct, c("barcode_rna", "cell_type"))
  
  # Read barcode mapping (ATAC = DNA, GEX = RNA)
  bmap <- fread(barcode_map_file, header = TRUE)
  setnames(bmap, c("barcode_dna", "barcode_rna"))
  
  # Merge: RNA barcode links cell type to DNA barcode
  merged <- merge(bmap, ct, by = "barcode_rna", all = FALSE)
  
  message(sprintf("  Mapped %d / %d RNA cell types to DNA barcodes (%d DNA barcodes in map)",
                  nrow(merged), nrow(ct), nrow(bmap)))
  
  result <- merged[, .(barcode = barcode_dna, cell_type)]
  return(result)
}

#' Groups cells by their cell type annotation, then performs hierarchical clustering within each group.
#'
#' @param mat Numeric matrix (bins x cells) used for clustering
#' @param cell_types data.table/data.frame with columns "barcode" and "cell_type"
#' @param dist_method Distance method (default: "manhattan")
#' @param hclust_method Agglomeration method (default: "ward.D2")
#' @return list with $ordered_cells (character vector), $cell_type_df (data.table
#'         with barcode, cell_type in plot order), $dendrograms (list of per-group
#'         dendrogram objects)
cluster_by_cell_type <- function(mat, cell_types,
                                 dist_method = "manhattan",
                                 hclust_method = "ward.D2") {
  # Ensure data.table
  ct <- as.data.table(cell_types)
  setnames(ct, c("barcode", "cell_type"))
  
  # Keep only cells present in the matrix
  available <- colnames(mat)
  ct <- ct[barcode %in% available]
  
  if (nrow(ct) == 0) {
    stop("No barcodes in cell_types match the profile column names.")
  }
  
  message(sprintf("  Matched %d / %d cells to cell type annotations",
                  nrow(ct), length(available)))
  
  # Cluster within each cell type group
  ordered_cells <- character(0)
  dendrograms <- list()
  
  for (ctype in sort(unique(ct$cell_type))) {
    group_cells <- ct[cell_type == ctype, barcode]
    
    if (length(group_cells) == 1) {
      # Single cell: no clustering needed
      ordered_cells <- c(ordered_cells, group_cells)
      dendrograms[[ctype]] <- NULL
    } else {
      # Hierarchical clustering within group
      sub_mat <- mat[, group_cells, drop = FALSE]
      sub_mat[is.na(sub_mat)] <- -1
      hc <- hclust(dist(t(sub_mat), method = dist_method), method = hclust_method)
      dhc <- as.dendrogram(hc)
      ordered_cells <- c(ordered_cells, labels(dhc))
      dendrograms[[ctype]] <- dhc
    }
  }
  
  # Build ordered cell_type data.frame
  cell_type_df <- ct[match(ordered_cells, barcode)]
  
  return(list(
    ordered_cells = ordered_cells,
    cell_type_df = cell_type_df,
    dendrograms = dendrograms
  ))
}

# -----------------------------------------------------------------------------
# Step 4: Data Extraction
# -----------------------------------------------------------------------------

# Extract genomic bins from ASCAT.sc results
#' @param res ASCAT.sc result object
#' @return data.table with chr, start, end, bin, start_cum, end_cum
extract_bins <- function(res) {
  message("Extracting genomic bins...")

  # Prefer nlSe (filtered by sc_excludeBadBins) over lSe (original full set).
  # nlSe uses $starts/$ends; lSe uses $start/$end.
  lSe_use    <- if (!is.null(res$nlSe)) res$nlSe else res$lSe
  use_plural <- !is.null(res$nlSe)

  bins_list <- lapply(names(lSe_use), function(chr) {
    data.table(
      chr   = chr,
      start = if (use_plural) lSe_use[[chr]]$starts else lSe_use[[chr]]$start,
      end   = if (use_plural) lSe_use[[chr]]$ends   else lSe_use[[chr]]$end
    )
  })
  
  # Combine all the chromosomes
  bins <- data.table::rbindlist(bins_list)
  
  # Add cumulative positions for genome-wide plotting
  # Give a sequential bin number across the genome
  bins[, bin := .I]
  
  # end_cum: cumulative position in bp (for x-axis)
  bins[, end_cum := cumsum((end - start) + 1)]
  
  # start_cum: starting cumulative position
  bins[, start_cum := c(1, end_cum[.I - 1] + 1)]
  
  message(sprintf("  Found %d bins across %d chromosomes", 
                  nrow(bins), length(unique(bins$chr))))
  
  return(bins)
}

# Make chromosome boundary table for plotting vertical lines
#' @param bins data.table from extract_bins()
#' @return data.table with chromosome boundaries
make_chr_bounds <- function(bins) {
  chr_bounds <- bins[, .(
    min = min(bin), 
    max = max(bin), 
    chrlen_bp = sum(end - start)
  ), by = chr]
  
  chr_bounds[, mid := round(min + (max - min) / 2, 0)]
  chr_bounds[, end_bp := cumsum(as.numeric(chrlen_bp))]
  chr_bounds[, start_bp := end_bp - chrlen_bp]
  chr_bounds[, mid_bp := round((chrlen_bp / 2) + start_bp, 0)]
  
  return(chr_bounds)
}

# Extract TOTAL copy number profiles (non-allele-specific)
#' Extracts both segmented CN and raw per-bin CN values
#' @param res ASCAT.sc result object
#' @return list with $profiles (data.table of CN per cell), $raw (raw CN per bin)
extract_total_profiles <- function(res) {
  message("Extracting total copy number profiles...")
  
  if (is.null(res$allProfiles)) {
    stop("Cannot find allProfiles in result object")
  }
  
  cell_names <- names(res$allProfiles)
  message(sprintf("  Processing %d cells...", length(cell_names)))
  
  # Get profiles
  profiles_list <- lapply(cell_names, function(cell) {
    dt <- as.data.table(res$allProfiles[[cell]])
    rep(dt$total_copy_number, as.numeric(dt$num.mark)) |> as.numeric()
  })
  
  profiles <- do.call(cbind, profiles_list) |> data.table()
  setnames(profiles, cell_names)
  
  # Get raw segments (per-bin smoothed values converted to CN scale)
  raw <- NULL
  if (!is.null(res$allTracks.processed) && !is.null(res$allSolutions)) {
    message("  Extracting raw per-bin CN values...")
    
    raw_list <- lapply(cell_names, function(cell) {
      dt <- rbindlist(res$allTracks.processed[[cell]]$lCTS)
      return(res$allSolutions[[cell]]$ploidy * 2^(dt$smoothed))
    })
    
    raw <- do.call(cbind, raw_list) |> data.table()
    setnames(raw, cell_names)
    message(sprintf("  Extracted raw data: %d bins x %d cells", nrow(raw), ncol(raw)))
  } else {
    message("  Warning: Could not find allTracks.processed or allSolutions")
    message("  Profile plots will not show raw data points")
  }
  
  message(sprintf("  Extracted %d cells x %d bins", ncol(profiles), nrow(profiles)))
  
  return(list(profiles = profiles, raw = raw))
}

# Extract ALLELE-SPECIFIC profiles (nMajor and nMinor)
#' Uses num.mark from res$allProfiles to expand segments to bins
#' @param res ASCAT.sc result object
#' @return list with $nMajor, $nMinor, $total (all data.tables)
extract_allele_profiles <- function(res, bins) {
  message("Extracting allele-specific profiles...")

  as_profiles <- if (!is.null(res$allProfiles_AS_smoothed)) res$allProfiles_AS_smoothed else res$allProfiles_AS
  if (is.null(as_profiles)) {
    stop("No allProfiles_AS or allProfiles_AS_smoothed found in result object.")
  }

  cell_names <- names(as_profiles)
  message(sprintf("  Processing %d cells...", length(cell_names)))
  
  n_bins <- nrow(bins)
  
  # allProfiles_AS_smoothed entries are data frames directly (nprof.fixed columns
  # + allele1Inferred/allele2Inferred/AS_mode_ON); allProfiles_AS entries are
  # lists with a $nprof.fixed element. Detect which structure we have.
  is_smoothed <- !is.null(res$allProfiles_AS_smoothed)

  profiles_list <- lapply(cell_names, function(cell) {
    cell_data <- as_profiles[[cell]]

    # Extract the profile data frame depending on structure
    if (is_smoothed) {
      # smoothed: cell_data IS the data frame
      if (!is.data.frame(cell_data) && !is.matrix(cell_data)) {
        return(list(nMajor = rep(NA_real_, n_bins),
                    nMinor = rep(NA_real_, n_bins),
                    total  = rep(NA_real_, n_bins)))
      }
      prof <- as.data.table(cell_data)
    } else {
      # unsmoothed: cell_data is a list with $nprof.fixed
      if (!is.list(cell_data) || is.null(cell_data$nprof.fixed)) {
        return(list(nMajor = rep(NA_real_, n_bins),
                    nMinor = rep(NA_real_, n_bins),
                    total  = rep(NA_real_, n_bins)))
      }
      prof <- as.data.table(cell_data$nprof.fixed)
    }

    # num.mark comes from the total (segmented) profiles — same row order
    total_prof <- as.data.table(res$allProfiles[[cell]])
    num_mark <- as.numeric(total_prof$num.mark)

    # If nA is NA, set to 0; if nB is NA, assign total_copy_number
    prof[is.na(nA), nA := 0]
    prof[is.na(nB), nB := total_copy_number]

    # Expand segments to bins using num.mark
    nA_expanded <- rep(as.numeric(prof$nA), num_mark)
    nB_expanded <- rep(as.numeric(prof$nB), num_mark)

    # nMajor = max(nA, nB), nMinor = min(nA, nB)
    nMajor <- pmax(nA_expanded, nB_expanded, na.rm = FALSE)
    nMinor <- pmin(nA_expanded, nB_expanded, na.rm = FALSE)
    total <- nMajor + nMinor

    list(nMajor = nMajor, nMinor = nMinor, total = total)
  })
  
  names(profiles_list) <- cell_names
  
  # Convert to data.tables (each column = one cell)
  nMajor <- as.data.table(do.call(cbind, lapply(profiles_list, `[[`, "nMajor")))
  nMinor <- as.data.table(do.call(cbind, lapply(profiles_list, `[[`, "nMinor")))
  total <- as.data.table(do.call(cbind, lapply(profiles_list, `[[`, "total")))
  
  setnames(nMajor, cell_names)
  setnames(nMinor, cell_names)
  setnames(total, cell_names)

  # Drop cells with no allele-specific call at all (every bin NA in both
  # alleles). These come from res$allProfiles_AS entries that were NULL or
  # lacked nprof.fixed -- ASCAT.sc didn't produce an AS decomposition for
  # them, so they would render as blank profiles and as a misleading flat
  # band in the heatmap.
  fully_na <- sapply(nMajor, function(col) all(is.na(col))) &
              sapply(nMinor, function(col) all(is.na(col)))
  dropped <- cell_names[fully_na]
  if (length(dropped) > 0) {
    message(sprintf("  Dropping %d cells with no allele-specific call (e.g. %s)",
                    length(dropped),
                    paste(head(dropped, 3), collapse = ", ")))
    keep <- !fully_na
    nMajor <- nMajor[, keep, with = FALSE]
    nMinor <- nMinor[, keep, with = FALSE]
    total  <- total[, keep, with = FALSE]
    cell_names <- cell_names[keep]
  }

  # Report NA counts across remaining cells (partial NA bins)
  na_per_cell <- sapply(nMajor, function(col) sum(is.na(col)))
  if (any(na_per_cell > 0)) {
    message(sprintf("  Warning: %d/%d cells have NA bins (max: %d bins = %.1f%%)",
                    sum(na_per_cell > 0), length(na_per_cell),
                    max(na_per_cell), 100 * max(na_per_cell) / n_bins))
  }

  message(sprintf("  Extracted %d cells x %d bins", ncol(nMajor), nrow(nMajor)))

  return(list(nMajor = nMajor, nMinor = nMinor, total = total,
              dropped_cells = dropped))
}

extract_baf_ci <- function(res, cell_name, bins) {
  as_profiles <- if (!is.null(res$allProfiles_AS_smoothed)) res$allProfiles_AS_smoothed else res$allProfiles_AS
  cell_as <- as_profiles[[cell_name]]

  # smoothed: cell_as IS the data frame; unsmoothed: access $nprof.fixed
  seg_prof <- as.data.table(if (is.data.frame(cell_as) || is.matrix(cell_as)) cell_as else cell_as$nprof.fixed)
  cell_total <- as.data.table(res$allProfiles[[cell_name]])
  num_mark <- as.numeric(cell_total$num.mark)
  
  total_cn <- as.numeric(seg_prof$total_copy_number)
  baf_lower <- as.numeric(seg_prof$q05)
  baf_upper <- as.numeric(seg_prof$q95)
  
  nA_int <- as.numeric(seg_prof$nA)
  nB_int <- as.numeric(seg_prof$nB)
  nA_int[is.na(nA_int)] <- 0
  nB_int[is.na(nB_int)] <- total_cn[is.na(nB_int)]
  
  # Selecting the major allele:
  nMajor_int <- pmax(nA_int, nB_int)
  nMinor_int <- pmin(nA_int, nB_int)
  
  # BAF quantiles scaled to copy-number space
  # ribbon reflects actual BAF uncertainty 
  nMajor_upper <- baf_upper * total_cn
  nMajor_lower <- baf_lower * total_cn
  nMinor_upper <- total_cn - nMajor_lower
  nMinor_lower <- total_cn - nMajor_upper
  
  # Expand to per-bin (include integer call so the plot can anchor the
  # ribbon between the call and the quantile bound)
  data.table(
    bin          = seq_len(nrow(bins)),
    nMajor       = rep(nMajor_int, num_mark),
    nMinor       = rep(nMinor_int, num_mark),
    nMajor_lower = rep(nMajor_lower, num_mark),
    nMajor_upper = rep(nMajor_upper, num_mark),
    nMinor_lower = rep(nMinor_lower, num_mark),
    nMinor_upper = rep(nMinor_upper, num_mark)
  )
}

# -----------------------------------------------------------------------------
# Step 5: Plotting Functions - Non-Allele-Specific
# -----------------------------------------------------------------------------

#' Create heatmap for TOTAL copy number (non-allele-specific)
#' @param profiles data.table of copy numbers (bins x cells)
#' @param bins data.table with genomic coordinates
#' @param dendrogram Show hierarchical clustering dendrogram
#' @param order Optional vector of cell names to order by
#' @param linesize Line width for heatmap
#' @param show_rownames Show cell names on y-axis
#' @param dist_method Distance method passed to \code{dist()} (default: "euclidean"; use "manhattan" for manhattan distance)
#' @param hclust_method Agglomeration method passed to \code{hclust()} (default: "average"; use "ward.D2" for Ward linkage)
#' @return ggplot object
# Function to plot genomewide copy number profile heatmaps

plotHeatmap = function(profiles, bins, order, dendrogram = TRUE,
                       linesize = .5, rownames = FALSE, annotation = NULL,
                       cell_types = NULL,
                       dist_method = "manhattan", hclust_method = "ward.D2") {
  # Check that order is not specified while dendrogram is also requested
  if (!missing(order)) dendrogram = FALSE
  
  # Make sure they are data.tables
  # Copy to avoid modifying caller's data by reference
  bins <- bins
  profiles <- copy(profiles)
  setDT(bins)
  setDT(profiles)
  
  # Set cowplot theme
  theme_set(theme_cowplot())
  # Get cumulative locations
  bins[, bin := seq_along(chr)]
  bins[, end_cum := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[seq_along(end_cum) - 1] + 1)]
  
  # Make chr_bounds
  chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end - start)), by = chr]
  chr_bounds = chr_bounds %>%
    mutate(mid = round(min + (max - min) / 2, 0),
           end_bp = cumsum(as.numeric(chrlen_bp)),
           start_bp = end_bp - chrlen_bp,
           mid_bp = round((chrlen_bp / 2) + start_bp, 0))
  
  # Colors
  colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
             "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
  
  names(colors) = c(as.character(0:10), "10+")
  
  dt = data.table(cbind(bins, profiles))
  
  # Set theme depending on rownames
  if (rownames) {
    custom_theme = theme(axis.ticks.y = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title = element_blank())
  } else {
    custom_theme = theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title = element_blank())
  }
  
  
  if (dendrogram) {
    # Clustering (use named columns instead of hardcoded index)
    bin_cols = c("chr", "start", "end", "bin", "start_cum", "end_cum")
    cell_cols = setdiff(names(dt), bin_cols)
    
    if (!is.null(cell_types)) {
      # --- Grouped clustering by cell type ---
      message("  Clustering within cell type groups...")
      mat <- as.matrix(dt[, ..cell_cols])
      ct_result <- cluster_by_cell_type(mat, cell_types,
                                        dist_method = dist_method,
                                        hclust_method = hclust_method)
      ordered_labels <- ct_result$ordered_cells
      
      # Build per-group dendrogram segments with y-offset so groups stack
      all_segments <- list()
      cumulative_offset <- 0
      
      for (ctype in sort(unique(ct_result$cell_type_df$cell_type))) {
        group_cells <- ct_result$cell_type_df[cell_type == ctype, barcode]
        n_group <- length(group_cells)
        
        if (n_group >= 2 && !is.null(ct_result$dendrograms[[ctype]])) {
          ddata_grp <- dendro_data(ct_result$dendrograms[[ctype]], type = "rectangle")
          seg_grp <- ggdendro::segment(ddata_grp)
          seg_grp$x    <- seg_grp$x    + cumulative_offset
          seg_grp$xend <- seg_grp$xend + cumulative_offset
          seg_grp$group <- ctype
          all_segments <- c(all_segments, list(seg_grp))
        }
        cumulative_offset <- cumulative_offset + n_group
      }
      
      n_samples <- length(ordered_labels)
      
      if (length(all_segments) > 0) {
        seg_df <- do.call(rbind, all_segments)
        dendro <- ggplot(seg_df) +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
          coord_flip() +
          scale_y_reverse(expand = c(0, 0)) +
          scale_x_continuous(limits = c(0.5, n_samples + 0.5), expand = c(0, 0)) +
          theme_dendro()
      } else {
        dendro <- ggplot() + theme_void()
      }
      
      # Build annotation from cell_types for the annotation bar
      annot_from_ct <- data.table(
        sample = ct_result$cell_type_df$barcode,
        variable = "Cell Type",
        value = ct_result$cell_type_df$cell_type
      )
      
    } else {
      # --- Standard global clustering ---
      hc = hclust(dist(t(dt[, ..cell_cols]), method = dist_method), method = hclust_method)
      dhc = as.dendrogram(hc)
      
      # Rectangular lines
      ddata = dendro_data(dhc, type = "rectangle")
      
      n_samples <- length(ddata$labels$label)
      ordered_labels <- ddata$labels$label
      dendro = ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 0)) +
        scale_x_continuous(limits = c(0.5, n_samples + 0.5), expand = c(0, 0)) +
        theme_dendro()
      
      annot_from_ct <- NULL
    }
    
    # Prepare for heatmap
    dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
    dt_melt[, value := as.character(value)]
    dt_melt[as.numeric(value) > 10, value := "10+"]
    dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
    
    # Add one invisible (NA-coordinate) row per absent CN level so that
    # geom_linerange always draws a legend key with the correct colour.
    .all_cn <- c(as.character(0:10), "10+")
    .missing_cn <- setdiff(.all_cn, levels(droplevels(dt_melt$value)))
    if (length(.missing_cn) > 0) {
      .var_levels <- levels(dt_melt$variable)
      dt_melt <- rbindlist(list(
        dt_melt,
        data.table(chr = dt_melt$chr[1], start = dt_melt$start[1],
                   end = dt_melt$end[1],  bin   = dt_melt$bin[1],
                   start_cum = NA_real_,  end_cum = NA_real_,
                   variable  = .var_levels[1], value = .missing_cn)
      ), use.names = TRUE, fill = TRUE)
      dt_melt[, value    := factor(as.character(value),    levels = .all_cn)]
      dt_melt[, variable := factor(as.character(variable), levels = .var_levels)]
    }
    
    # Set sample order
    dt_melt[, variable  := factor(variable, levels = ordered_labels)]
    dt_melt[, mid_cum   := (start_cum + end_cum) / 2]
    dt_melt[, bin_width := end_cum - start_cum]
    
    # Plot heatmap
    heatmap = ggplot(dt_melt) +
      geom_tile(aes(x = mid_cum, y = variable, fill = value, width = bin_width), height = 1) +
      scale_fill_manual(values = colors, drop = FALSE) +
      labs(fill = "Copy Number") +
      scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) +
      geom_vline(data = chr_bounds, aes(xintercept = end_bp), linetype = 1, linewidth = .8) +
      custom_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    # Use cell_types annotation if provided, otherwise fall back to annotation param
    use_annotation <- if (!is.null(annot_from_ct)) annot_from_ct else annotation
    
    if (!is.null(use_annotation)) {
      setnames(use_annotation, c("sample", "variable", "value"))
      
      # Plot annotation
      use_annotation[, sample := factor(sample, levels = ordered_labels)]
      annot_plt = ggplot(use_annotation, aes(x = 1, y = sample, fill = value)) +
        geom_bar(stat = "identity", width = 1) +
        scale_fill_npg() +
        theme_void() +
        theme(legend.position = "right")
      
      combined = dendro + annot_plt + heatmap + plot_layout(ncol = 3, widths = c(0.2, 0.05, 2), guides = "collect")
      
    } else {
      combined =cowplot::plot_grid(dendro, heatmap, ncol = 2, rel_widths = c(0.1, 1), align = "h", axis = "tb")
    }
    
    return(combined)
  }
  
  # Prepare for heatmap
  dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
  dt_melt[, value := as.character(value)]
  dt_melt[as.numeric(value) > 10, value := "10+"]
  dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
  
  # Add one invisible (NA-coordinate) row per absent CN level so that
  # geom_linerange always draws a legend key with the correct colour.
  .all_cn <- c(as.character(0:10), "10+")
  .missing_cn <- setdiff(.all_cn, levels(droplevels(dt_melt$value)))
  if (length(.missing_cn) > 0) {
    .var_levels <- levels(dt_melt$variable)
    dt_melt <- rbindlist(list(
      dt_melt,
      data.table(chr = dt_melt$chr[1], start = dt_melt$start[1],
                 end = dt_melt$end[1],  bin   = dt_melt$bin[1],
                 start_cum = NA_real_,  end_cum = NA_real_,
                 variable  = .var_levels[1], value = .missing_cn)
    ), use.names = TRUE, fill = TRUE)
    dt_melt[, value    := factor(as.character(value),    levels = .all_cn)]
    dt_melt[, variable := factor(as.character(variable), levels = .var_levels)]
  }
  
  if (!missing(order)) {
    # Set sample order
    dt_melt[, variable := factor(variable, levels = order)]
  } else {
    # Clustering (use named columns instead of hardcoded index)
    bin_cols = c("chr", "start", "end", "bin", "start_cum", "end_cum")
    cell_cols = setdiff(names(dt), bin_cols)
    hc = hclust(dist(t(dt[, ..cell_cols]), method = dist_method), method = hclust_method)
    dhc = as.dendrogram(hc)
    
    dt_melt[, variable := factor(variable, levels = labels(dhc))]
  }
  dt_melt[, mid_cum   := (start_cum + end_cum) / 2]
  dt_melt[, bin_width := end_cum - start_cum]
  
  # Plot heatmap
  heatmap = ggplot(dt_melt) +
    geom_tile(aes(x = mid_cum, y = variable, fill = value, width = bin_width), height = 1) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(fill = "Copy Number") +
    scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) +
    geom_vline(data = chr_bounds, aes(xintercept = end_bp), linetype = 1, linewidth = .8) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  return(heatmap)
}

# Create single-cell profile plot for TOTAL copy number
plotProfile = function(segments, raw, bins, sc = TRUE, linesize = 1, cell_name = NULL, ploidy = NULL,
                       sex = c("auto", "male", "female")) {
  sex <- match.arg(sex)
  
  # Set theme
  theme_set(theme_cowplot())
  
  # Process bins
  bins = bins
  bins[, bin := seq_along(chr)]
  bins[, end_cum := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[seq_along(end_cum)-1] + 1)]
  
  # Make chr_bounds
  chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end - start)), by = chr]
  chr_bounds = chr_bounds %>%
    mutate(mid = round(min + (max - min) / 2, 0),
           end_bp = cumsum(as.numeric(chrlen_bp)),
           start_bp = end_bp - chrlen_bp,
           mid_bp = round((chrlen_bp / 2) + start_bp, 0))
  
  if (sc) {
    # dt = cbind(bins, segments, raw * cn)
    dt = cbind(bins, segments, raw)
    setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "cn", "raw"))
    
    # ---- Sex chromosome handling -----------------------------------------
    sex_chr_lvls <- c("chrX", "chrY", "X", "Y", "23", "24")
    is_male <- sex == "male"
    if (sex == "auto") {
      chrY_cn <- dt[chr %in% c("chrY", "Y", "24"), cn]
      chrY_cn <- chrY_cn[!is.na(chrY_cn)]
      if (length(chrY_cn) > 0 && stats::median(chrY_cn) >= 1) is_male <- TRUE
    }
    
    # Cap CN at 6 for visual display – keeps a fixed 0-6 axis across all
    # profiles so they remain directly comparable. Anything > 5 is "5+".
    dt[, cn_display := pmin(cn, 6)]
    dt[, col_cn := pmin(cn, 6)]
    if (is_male) {
      message("  Male sample detected: shifting chrX/chrY colour baseline (+1)")
      dt[chr %in% sex_chr_lvls, col_cn := pmin(cn + 1L, 6L)]
    }
    dt[, col := ifelse(col_cn < 6, as.character(col_cn), "5+")]
    dt[, col := factor(col, levels = c(as.character(0:5), "5+"))]
    
    colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e", "#821010")
    names(colors) = c(as.character(0:5), "5+")
    
    # save plot
    plot = ggplot(dt, aes(x = bin)) +
      geom_point(aes(y = raw, color = col), size = 0.7) +
      geom_segment(aes(x = bin - 0.5, xend = bin + 0.5, y = cn_display, yend = cn_display), linewidth = linesize) +
      scale_color_manual(values = colors, drop = FALSE) +
      scale_y_continuous(labels = scales::comma_format(accuracy = 1),
                         breaks = c(0, 2, 4, 6),
                         limits = c(0, 6),
                         oob = scales::squish) +
      scale_x_continuous(expand = c(0, 0), breaks = chr_bounds$mid, labels = chr_bounds$chr) +
      labs(y = "Copy Number", x = "") +
      geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
            axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 16))
  } else {
    dt = cbind(bins, segments, raw)
    setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "segment", "raw"))
    
    # Assign colors
    dt[segment > log2(2.5 / 2), col := "Amplification"]
    dt[segment < log2(1.5 / 2), col := "Deletion"]
    dt[is.na(col), col := "Neutral"]
    
    colors = c(brewer.pal(3, "Set1")[1:2], "#c1c1c1")
    names(colors) = c("Amplification", "Deletion", "Neutral")
    
    # save plot
    plot = ggplot(dt, aes(x = bin)) +
      geom_point(aes(y = raw, color = col), size = 0.7) +
      geom_point(aes(y = segment), size = linesize) +
      scale_color_manual(values = colors, drop = FALSE) +
      scale_y_continuous(limits = c(-4, 4), labels = scales::comma_format(accuracy = 1), breaks = scales::pretty_breaks(6)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(y = "Log2 ratio", x = "") +
      geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
      geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", angle = 45, inherit.aes = FALSE) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  # Add title with ploidy
  if (!is.null(cell_name)) {
    title_text <- if (!is.null(ploidy)) {
      sprintf("%s (ploidy = %.2f)", cell_name, ploidy)
    } else {
      cell_name
    }
    plot <- plot + ggtitle(title_text) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
  }
  
  return(plot)
}

# -----------------------------------------------------------------------------
# Step 6: Plotting Functions - Allele-Specific
# -----------------------------------------------------------------------------

#' Build custom allele-specific legend grouped by total CN
#' Matches the aesthetic of ASCAT-style publications with gradient swatches
#' @return ggplot object (legend strip)
build_allele_legend <- function() {
  as_colors <- get_allele_specific_colors()
  
  # Define states grouped by total CN (excluding 0|0 which is homozyg deletion)
  legend_states <- list(
    "1" = data.table(state = "1|0", label = "(1, 0)"),
    "2" = data.table(state = c("1|1", "2|0"), label = c("(1, 1)", "(2, 0)")),
    "3" = data.table(state = c("2|1", "3|0"), label = c("(2, 1)", "(3, 0)")),
    "4" = data.table(state = c("2|2", "3|1", "4|0"), label = c("(2, 2)", "(3, 1)", "(4, 0)")),
    "5" = data.table(state = c("3|2", "4|1", "5|0"), label = c("(3, 2)", "(4, 1)", "(5, 0)")),
    "6" = data.table(state = c("3|3", "4|2", "5|1", "6|0"), label = c("(3, 3)", "(4, 2)", "(5, 1)", "(6, 0)")),
    "7" = data.table(state = c("4|3", "5|2", "6|1", "7|0"), label = c("(4, 3)", "(5, 2)", "(6, 1)", "(7, 0)"))
  )
  
  # Build data.table with x positions (tiles side-by-side, gaps between groups)
  tile_width <- 1
  gap <- 0.8
  rows <- list()
  current_x <- 0
  
  for (cn in names(legend_states)) {
    grp <- legend_states[[cn]]
    n <- nrow(grp)
    grp[, total_cn := cn]
    grp[, x := current_x + (seq_len(n) - 1) * tile_width + tile_width / 2]
    grp[, x_center := current_x + (n * tile_width) / 2]
    rows <- c(rows, list(grp))
    current_x <- current_x + n * tile_width + gap
  }
  
  legend_data <- rbindlist(rows)
  
  # Group header positions (one per total CN)
  group_headers <- legend_data[, .(x_mid = x_center[1],
                                   header = paste0("Total CN=", total_cn[1])),
                               by = total_cn]
  
  p <- ggplot(legend_data) +
    # Colored tiles
    geom_tile(aes(x = x, y = 0, fill = state),
              width = tile_width * 0.92, height = 0.7,
              color = "grey30", linewidth = 0.15) +
    scale_fill_manual(values = as_colors, guide = "none") +
    # State labels below tiles
    geom_text(aes(x = x, y = -0.65, label = label),
              size = 2.0, color = "grey20") +
    # Group headers above tiles
    geom_text(data = group_headers,
              aes(x = x_mid, y = 0.7, label = header),
              size = 2.5, fontface = "bold") +
    # "Copy number (CN)" label on the left
    annotate("text", x = -1.2, y = 0, label = "Copy number\n(CN)",
             size = 2.5, fontface = "bold", hjust = 1, lineheight = 0.9) +
    coord_cartesian(clip = "off",
                    xlim = c(-2.5, current_x),
                    ylim = c(-1.2, 1.1)) +
    theme_void() +
    theme(plot.margin = margin(2, 10, 2, 35, unit = "pt"))
  
  return(p)
}

#' Create heatmap for ALLELE-SPECIFIC copy number
#' @param nMajor data.table of major allele CN (bins x cells)
#' @param nMinor data.table of minor allele CN (bins x cells)
#' @param bins data.table with genomic coordinates
#' @param dendrogram Show hierarchical clustering dendrogram
#' @param order Optional vector of cell names to order by
#' @param linesize Line width for heatmap
#' @param show_rownames Show cell names on y-axis
#' @param title Optional plot title
#' @return ggplot object

plot_allele_heatmap <- function(nMajor, nMinor, bins,
                                dendrogram = TRUE,
                                order = NULL,
                                show_rownames = FALSE,
                                cell_types = NULL,
                                title = NULL) {
  # Copy to avoid modifying caller's data by reference
  bins <- bins
  nMajor <- copy(nMajor)
  nMinor <- copy(nMinor)
  setDT(bins)
  setDT(nMajor)
  setDT(nMinor)
  
  # === FIX FOR DUPLICATE CELL NAMES ===
  orig_names <- colnames(nMajor)
  if (any(duplicated(orig_names))) {
    warning(sprintf("Found %d duplicate cell names, making unique", sum(duplicated(orig_names))))
    new_names <- make.unique(orig_names)
    setnames(nMajor, new_names)
    setnames(nMinor, new_names)
  }
  # === END FIX ===
  
  # Get colors and chromosome bounds
  as_colors <- get_allele_specific_colors()
  # Ensure "NA" is included in state levels
  state_levels <- names(as_colors)
  chr_bounds <- make_chr_bounds(bins)
  
  # Initialize
  dendro_plot <- NULL
  ordered_samples <- colnames(nMajor)
  
  # Handle ordering
  if (!is.null(order)) {
    ordered_samples <- order[order %in% colnames(nMajor)]
    dendrogram <- FALSE
  }
  
  annot_from_ct <- NULL
  
  if (dendrogram) {
    message("  Performing hierarchical clustering on allele-specific states...")
    
    # Use Manhattan distance on actual nMajor/nMinor values
    maj_mat <- as.matrix(nMajor)
    min_mat <- as.matrix(nMinor)
    combined <- rbind(maj_mat, min_mat)
    
    # Replace NA with -1 for clustering (NA becomes a distinct state)
    combined[is.na(combined)] <- -1
    
    if (!is.null(cell_types)) {
      # --- Grouped clustering by cell type ---
      message("  Clustering within cell type groups...")
      ct_result <- cluster_by_cell_type(combined, cell_types,
                                        dist_method = "manhattan",
                                        hclust_method = "ward.D2")
      ordered_samples <- ct_result$ordered_cells
      
      # Build per-group dendrogram segments with y-offset
      all_segments <- list()
      cumulative_offset <- 0
      
      for (ctype in sort(unique(ct_result$cell_type_df$cell_type))) {
        group_cells <- ct_result$cell_type_df[cell_type == ctype, barcode]
        n_group <- length(group_cells)
        
        if (n_group >= 2 && !is.null(ct_result$dendrograms[[ctype]])) {
          ddata_grp <- dendro_data(ct_result$dendrograms[[ctype]], type = "rectangle")
          seg_grp <- ggdendro::segment(ddata_grp)
          seg_grp$x    <- seg_grp$x    + cumulative_offset
          seg_grp$xend <- seg_grp$xend + cumulative_offset
          seg_grp$group <- ctype
          all_segments <- c(all_segments, list(seg_grp))
        }
        cumulative_offset <- cumulative_offset + n_group
      }
      
      n_samples <- length(ordered_samples)
      
      if (length(all_segments) > 0) {
        seg_df <- do.call(rbind, all_segments)
        dendro_plot <- ggplot(seg_df) +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
          coord_flip() +
          scale_y_reverse(expand = c(0, 0)) +
          scale_x_continuous(limits = c(0.5, n_samples + 0.5), expand = c(0, 0)) +
          theme_dendro()
      } else {
        dendro_plot <- ggplot() + theme_void()
      }
      
      # Build annotation from cell_types
      annot_from_ct <- data.table(
        sample = ct_result$cell_type_df$barcode,
        variable = "Cell Type",
        value = ct_result$cell_type_df$cell_type
      )
      
    } else {
      # --- Standard global clustering ---
      hc <- hclust(dist(t(combined), method = "manhattan"), method = "ward.D2")
      dhc <- as.dendrogram(hc)
      ordered_samples <- labels(dhc)
      
      # Create dendrogram plot - matching plotHeatmap approach exactly
      ddata <- dendro_data(dhc, type = "rectangle")
      n_samples <- length(ordered_samples)
      
      dendro_plot <- ggplot(segment(ddata)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 0)) +
        scale_x_continuous(limits = c(0.5, n_samples + 0.5), expand = c(0, 0)) +
        theme_dendro()
    }
  }
  
  # Subset and order profiles
  nMajor_ordered <- nMajor[, ..ordered_samples]
  nMinor_ordered <- nMinor[, ..ordered_samples]
  
  # Create allele-state matrix
  message("  Encoding allele-specific states...")
  
  # Store the numeric coordinates separately to avoid data.table type coercion issues
  coord_cols <- bins[, .(chr, start, end, bin, start_cum, end_cum)]
  
  # Create a separate data.table for states only
  state_list <- lapply(ordered_samples, function(cell) {
    encode_allele_state(nMajor_ordered[[cell]], nMinor_ordered[[cell]])
  })
  names(state_list) <- ordered_samples
  state_dt <- as.data.table(state_list)
  
  # Combine coordinates and states
  allele_states <- cbind(coord_cols, state_dt)
  
  # Melt to long format
  dt_melt <- melt(allele_states,
                  id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"),
                  variable.name = "cell",
                  value.name = "state")
  
  # Ensure numeric columns remain numeric after melt
  dt_melt[, start_cum := as.numeric(start_cum)]
  dt_melt[, end_cum := as.numeric(end_cum)]
  dt_melt[, mid_cum := (start_cum + end_cum) / 2]
  dt_melt[, bin_width := end_cum - start_cum]
  
  # Set cell as factor with ordered_samples levels (matching plotHeatmap approach)
  dt_melt[, cell := factor(as.character(cell), levels = ordered_samples)]
  
  # Factor with proper levels for state
  dt_melt[, state := factor(state, levels = state_levels)]
  
  # Theme
  custom_theme <- if (show_rownames) {
    theme(axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 6))
  } else {
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank())
  }
  
  heatmap <- ggplot(dt_melt) +
    geom_tile(aes(x = mid_cum, y = cell, fill = state, width = bin_width), height = 1) +
    scale_fill_manual(values = as_colors, limits = names(as_colors), breaks = names(as_colors),
                      na.value = "#BEBEBE", guide = "none") +
    scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_vline(data = chr_bounds, aes(xintercept = end_bp),
               linetype = 1, linewidth = 0.3) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Add title
  if (!is.null(title)) {
    heatmap <- heatmap + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  }
  
  # Build custom allele-specific legend
  legend_plot <- build_allele_legend()
  
  # Combine with dendrogram and optional annotation bar
  if (!is.null(dendro_plot) && !is.null(annot_from_ct)) {
    # With dendrogram + cell type annotation bar
    annot_from_ct[, sample := factor(sample, levels = ordered_samples)]
    annot_plt <- ggplot(annot_from_ct, aes(x = 1, y = sample, fill = value)) +
      geom_bar(stat = "identity", width = 1) +
      scale_fill_npg() +
      theme_void() +
      theme(legend.position = "right")
    
    main_plot <- dendro_plot + annot_plt + heatmap +
      plot_layout(ncol = 3, widths = c(0.2, 0.05, 2), guides = "collect")
  } else if (!is.null(dendro_plot)) {
    main_plot <- cowplot::plot_grid(dendro_plot, heatmap, ncol = 2,
                                    rel_widths = c(0.1, 1), align = "h", axis = "tb")
  } else {
    main_plot <- heatmap
  }
  
  # Stack main plot and custom legend vertically
  final <- main_plot / legend_plot + plot_layout(heights = c(1, 0.06))
  
  return(final)
}

#' Create single-cell profile plot for ALLELE-SPECIFIC data
#' @param nMajor Vector of major allele CN
#' @param nMinor Vector of minor allele CN
#' @param bins data.table with genomic coordinates
#' @param cell_name Name for the title
#' @param bar_height Height of each allele bar (default 0.15)
#' @return ggplot object
plot_allele_profile <- function(nMajor, nMinor, bins, cell_name = NULL, ploidy = NULL, segments_ci = NULL, bar_height = 0.05, ci_alpha = 0.2) {
  
  # Ensure major >= minor consistently, same as extract_allele_profiles
  nMaj_orig <- nMajor
  nMin_orig <- nMinor
  nMajor <- pmax(nMaj_orig, nMin_orig, na.rm = FALSE)
  nMinor <- pmin(nMaj_orig, nMin_orig, na.rm = FALSE)
  
  # Set theme to match plotProfile
  theme_set(theme_cowplot())
  
  setDT(bins)
  bins[, bin := seq_along(chr)]
  bins[, end_cum := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[seq_along(end_cum) - 1] + 1)]
  
  # Make chr_bounds - same method as plotProfile
  chr_bounds <- bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end - start)), by = chr]
  chr_bounds <- chr_bounds %>%
    mutate(mid = round(min + (max - min) / 2, 0),
           end_bp = cumsum(as.numeric(chrlen_bp)),
           start_bp = end_bp - chrlen_bp,
           mid_bp = round((chrlen_bp / 2) + start_bp, 0))
  
  # Create segment data (consecutive bins with same CN get merged)
  dt <- data.table(
    bin = bins$bin,
    chr = bins$chr,
    nMajor = nMajor,
    nMinor = nMinor
  )
  
  # Create segments by detecting CN changes
  dt[, seg_id := rleid(chr, nMajor, nMinor)]
  
  # Summarize segments
  segments <- dt[, .(
    start_bin = min(bin),
    end_bin = max(bin),
    nMajor = nMajor[1],
    nMinor = nMinor[1],
    chr = chr[1]
  ), by = seg_id]
  
  # Calculate rectangle positions
  # Major: sits on top of the CN line
  # Minor: sits below the CN line
  segments[, maj_ymin := nMajor]
  segments[, maj_ymax := nMajor + bar_height]
  segments[, min_ymin := nMinor - bar_height]
  segments[, min_ymax := nMinor]
  
  # Plot - ASCAT style with rectangles, matching plotProfile styling
  p <- ggplot(segments) +
    # Horizontal gridlines at integer CN values (behind the data)
    geom_hline(yintercept = 0:4, color = "grey80", linewidth = 0.3)
  if (!is.null(segments_ci)) {
    ci_dt <- copy(segments_ci)[!is.na(nMajor_lower)]
    
    # Fill from the integer call out to each quantile bound so the ribbon
    # reads as "call sits here, BAF uncertainty extends this far"
    ci_dt[, `:=`(
      maj_ribbon_lo = pmin(nMajor, nMajor_lower),
      maj_ribbon_hi = pmax(nMajor, nMajor_upper),
      min_ribbon_lo = pmin(nMinor, nMinor_lower),
      min_ribbon_hi = pmax(nMinor, nMinor_upper)
    )]

    p <- p +
      geom_rect(data = ci_dt,
                aes(xmin = bin - 0.5, xmax = bin + 0.5,
                    ymin = maj_ribbon_lo, ymax = maj_ribbon_hi),
                fill = "#E03546", alpha = ci_alpha, color = NA) +
      geom_rect(data = ci_dt,
                aes(xmin = bin - 0.5, xmax = bin + 0.5,
                    ymin = min_ribbon_lo, ymax = min_ribbon_hi),
                fill = "#1b38ae", alpha = ci_alpha, color = NA)
  }
  p <- p +
    # Major allele (red) - rectangle above CN line
    geom_rect(aes(xmin = start_bin - 0.5, xmax = end_bin + 0.5,
                  ymin = maj_ymin, ymax = maj_ymax,
                  fill = "Major"),
              color = NA) +
    # Minor allele (blue) - rectangle below CN line
    geom_rect(aes(xmin = start_bin - 0.5, xmax = end_bin + 0.5,
                  ymin = min_ymin, ymax = min_ymax,
                  fill = "Minor"),
              color = NA) +
    # Color scale with legend
    scale_fill_manual(
      values = c("Major" = "#E03546", "Minor" = "#1b38ae"),
      name = "Allele"
    ) +
    # Chromosome boundaries - dashed lines like plotProfile
    geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
    # Axis settings - allele-specific range (0-4 per allele, integer breaks
    # so a diploid balanced state shows both alleles at CN = 1)
    scale_y_continuous(
      labels = scales::comma_format(accuracy = 1),
      breaks = 0:4,
      limits = c(-bar_height, 4),
      oob = scales::squish,
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = chr_bounds$mid,
      labels = chr_bounds$chr
    ) +
    labs(y = "Copy Number", x = "") +
    # Theme matching plotProfile
    theme(
      legend.position = "top",
      legend.justification = "left",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 16)
    )
  
  # Add title with ploidy
  if (!is.null(cell_name)) {
    title_text <- if (!is.null(ploidy)) {
      sprintf("%s (ploidy = %.2f)", cell_name, ploidy)
    } else {
      cell_name
    }
    p <- p + ggtitle(title_text) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
  }
  
  return(p)
}

# -----------------------------------------------------------------------------
# STEP 6b: Auto-Detection Helper
# -----------------------------------------------------------------------------

#' Detect if ASCAT.sc result has allele-specific data
#' @param res ASCAT.sc result object
#' @return Logical: TRUE if allele-specific, FALSE otherwise
detect_allele_specific <- function(res) {
  if (!is.null(res$allProfiles_AS_smoothed) && length(res$allProfiles_AS_smoothed) > 0) {
    return(TRUE)
  }
  if (!is.null(res$allProfiles_AS) && length(res$allProfiles_AS) > 0) {
    return(TRUE)
  }
  
  first_cell <- names(res$allProfiles)[1]
  first_profile <- res$allProfiles[[first_cell]]
  
  if (is.data.frame(first_profile)) {
    return(all(c("nMajor", "nMinor") %in% colnames(first_profile)))
  }
  
  return(FALSE)
}


#' Print summary of ASCAT.sc result structure
#' @param res ASCAT.sc result object
summarize_ascat_result <- function(res) {
  cat("=== ASCAT.sc Result Summary ===\n")
  cat(sprintf("Elements: %s\n", paste(names(res), collapse = ", ")))
  
  if (!is.null(res$allProfiles)) {
    n_cells <- length(res$allProfiles)
    first_cell <- names(res$allProfiles)[1]
    first_profile <- as.data.frame(res$allProfiles[[first_cell]])
    
    cat(sprintf("Number of cells: %d\n", n_cells))
    cat(sprintf("Profile columns: %s\n", paste(colnames(first_profile), collapse = ", ")))
    
    is_allele_specific <- all(c("nMajor", "nMinor") %in% colnames(first_profile))
    cat(sprintf("Data type: %s\n", 
                ifelse(is_allele_specific, "ALLELE-SPECIFIC", "TOTAL COPY NUMBER ONLY")))
    
    # Calculate total bins
    total_bins <- sum(first_profile$num.mark)
    cat(sprintf("Total bins: %d\n", total_bins))
  }
  
  if (!is.null(res$build)) cat(sprintf("Genome build: %s\n", res$build))
  if (!is.null(res$binsize)) cat(sprintf("Bin size: %s\n", res$binsize))
  
  cat("===============================\n")
}

# -----------------------------------------------------------------------------
# STEP 7: Save Function - Outputs PNG and/or PDF
# -----------------------------------------------------------------------------

#' Save plot to PNG and/or PDF
#' @param plot ggplot object
#' @param output_dir Directory for output files
#' @param output_prefix Filename prefix (without extension)
#' @param output_format Vector of formats: "png", "pdf", or both
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param dpi Resolution for PNG
save_plot <- function(plot, 
                      output_dir = ".", 
                      output_prefix = "ASCAT_plot",
                      output_format = c("png", "pdf"),
                      width = 20,
                      height = 10,
                      dpi = 300) {
  
  # Create directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  base_path <- file.path(output_dir, output_prefix)
  
  # Save PNG
  if ("png" %in% output_format) {
    png_file <- paste0(base_path, ".png")
    
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggsave(png_file, plot, width = width, height = height, 
             dpi = dpi, device = ragg::agg_png, bg = "white")
    } else {
      # Fallback to standard PNG
      ggsave(png_file, plot, width = width, height = height, 
             dpi = dpi, device = "png", bg = "white")
    }
    message(sprintf("  Saved: %s", png_file))
  }
  
  # Save PDF
  if ("pdf" %in% output_format) {
    pdf_file <- paste0(base_path, ".pdf")
    ggsave(pdf_file, plot, width = width, height = height, device = "pdf")
    message(sprintf("  Saved: %s", pdf_file))
  }
}

# =============================================================================
# STEP 8: MAIN FUNCTION - The User-Facing Interface
# =============================================================================

#' Plot ASCAT.sc results
#' 
#' Main function for visualizing ASCAT.sc copy number results.
#' Supports both allele-specific and total copy number visualization.
#' 
#' @param res ASCAT.sc result object (from readRDS)
#' @param allele_specific Logical or "auto". If TRUE, use allele-specific coloring.
#'        If FALSE, use total copy number coloring. If "auto", detect from data.
#' @param plot_type Character. "all" (default) for both heatmap and all profiles,
#'        "heatmap" for genome-wide heatmap only, "profile" for profile plots only.
#' @param cell_id Character vector. For plot_type="profile", which cells to plot.
#'        If NULL and plot_type="all", plots all cells. If NULL and plot_type="profile",
#'        plots only the first cell.
#' @param output_dir Directory for saving output files.
#' @param output_prefix Filename prefix for heatmap (without extension).
#' @param output_format Output format: "png" (default).
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi Resolution for PNG output.
#' @param dendrogram Logical. Show clustering dendrogram in heatmap.
#' @param order Optional vector of cell names for custom ordering.
#' @param linesize Line width for heatmap / segment lines.
#' @param show_rownames Logical. Show cell names on y-axis.
#' @param cell_types Optional data.frame/data.table with columns "barcode" and "cell_type".
#'        When provided, cells are grouped by cell type and hierarchically clustered
#'        (manhattan distance, ward.D2) within each group. An annotation color bar is added.
#' @param title Optional plot title.
#' @param sex Sample sex: "auto" (default), "male", or "female". Only affects the
#'        TOTAL-CN profile plots: in male samples chrX/chrY are coloured with a
#'        baseline of CN=1 (so CN=1 reads as neutral grey instead of loss). Auto
#'        detection flags a profile as male when median chrY CN is >= 1.
#' @param save Logical. If TRUE, save to files. If FALSE, just return plot.
#' 
#' @return ggplot object (invisibly if save=TRUE)
#' 
#' @examples
#' # Load data
#' res <- readRDS("ASCAT_results.rds")
#' 
#' # Generate all plots (heatmap + all cell profiles)
#' ascatsc_plot(res)
#' 
#' # Only heatmap
#' ascatsc_plot(res, plot_type = "heatmap")
#' 
#' # Only profiles for specific cells
#' ascatsc_plot(res, plot_type = "profile", cell_id = c("AGCTAGTAGGCCTTGT", "AAACAAGCAAAGTCAT"))
#'
#' # Heatmap grouped by cell type with within-group clustering
#' # First map RNA cell types to DNA barcodes
#' ct <- map_cell_types_to_dna(
#'   cell_type_file = "final_cell_type.txt",
#'   barcode_map_file = "barcodes_atac_gex.csv"
#' )
#' ascatsc_plot(res, plot_type = "heatmap", cell_types = ct)
#'
ascatsc_plot <- function(res,
                         allele_specific = "auto",
                         plot_type = c("all", "heatmap", "profile"),
                         cell_id = NULL,
                         output_dir = ".",
                         output_prefix = "heatmap",
                         output_format = c("png"),
                         width = 20,
                         height = 10,
                         dpi = 300,
                         dendrogram = TRUE,
                         order = NULL,
                         linesize = 0.5,
                         show_ci = FALSE,
                         show_rownames = FALSE,
                         cell_types = NULL,
                         title = NULL,
                         sex = c("auto", "male", "female"),
                         save = TRUE) {
  
  sex <- match.arg(sex)
  
  # Match arguments
  
  plot_type <- match.arg(plot_type)
  
  # Handle auto-detection of allele_specific
  if (identical(allele_specific, "auto")) {
    allele_specific <- detect_allele_specific(res)
    message(sprintf("Auto-detected data type: %s", 
                    ifelse(allele_specific, "Allele-Specific", "Total Copy Number")))
  }
  
  message("=== ASCAT.sc Plotting ===")
  message(sprintf("Mode: %s", ifelse(allele_specific, "Allele-Specific", "Total Copy Number")))
  message(sprintf("Plot type: %s", plot_type))
  
  # -------------------------------------------------------------------------
  # STEP A: Extract bins (common to both modes)
  # -------------------------------------------------------------------------
  bins <- extract_bins(res)
  
  # -------------------------------------------------------------------------
  # STEP B: Create output directories
  # -------------------------------------------------------------------------
  plots_dir <- file.path(output_dir, "plots")
  profiles_dir <- file.path(plots_dir, "profiles")
  
  if (save) {
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    if (!dir.exists(profiles_dir)) dir.create(profiles_dir, recursive = TRUE)
  }
  
  # -------------------------------------------------------------------------
  # STEP C: Extract profiles and create plots
  # -------------------------------------------------------------------------
  if (allele_specific) {
    # Allele-specific mode
    profiles <- extract_allele_profiles(res, bins)
    cell_names <- colnames(profiles$nMajor)
    
    # Subset to requested cells if cell_id is provided
    if (!is.null(cell_id)) {
      matching <- cell_id[cell_id %in% cell_names]
      if (length(matching) == 0) {
        warning("No cell_id matched. Available names look like: ",
                paste(head(cell_names, 3), collapse = ", "))
      } else {
        message(sprintf("  Subsetting to %d / %d requested cells",
                        length(matching), length(cell_id)))
        profiles$nMajor <- profiles$nMajor[, ..matching]
        profiles$nMinor <- profiles$nMinor[, ..matching]
        profiles$total  <- profiles$total[, ..matching]
        cell_names <- matching
      }
    }
    
    # Create heatmap if requested
    if (plot_type %in% c("all", "heatmap")) {
      message("Creating allele-specific heatmap...")
      
      heatmap_plot <- plot_allele_heatmap(
        nMajor = profiles$nMajor,
        nMinor = profiles$nMinor,
        bins = bins,
        dendrogram = dendrogram,
        order = order,
        show_rownames = show_rownames,
        cell_types = cell_types,
        title = title
      )
      
      if (save) {
        # heatmap_file <- file.path(plots_dir, paste0(output_prefix, ".png"))
        ext <- output_format[1]
        heatmap_file <- file.path(plots_dir, paste0(output_prefix, ".", ext))
        ggsave(heatmap_file, heatmap_plot, width = width, height = height, dpi = dpi, bg = "white")
        message(sprintf("  Saved: %s", heatmap_file))
      }
    }
    
    # Create profile plots if requested
    if (plot_type %in% c("all", "profile")) {
      # Determine which cells to plot
      if (!is.null(cell_id)) {
        cells_to_plot <- cell_id
      } else if (plot_type == "all") {
        cells_to_plot <- cell_names
      } else {
        cells_to_plot <- cell_names[1]
        message(sprintf("No cell_id specified, using first cell: %s", cells_to_plot))
      }
      
      message(sprintf("Creating profile plots for %d cells...", length(cells_to_plot)))
      
      for (cell in cells_to_plot) {
        if (!cell %in% cell_names) {
          warning(sprintf("Cell '%s' not found, skipping", cell))
          next
        }
        
        ci_data <- NULL
        if (show_ci) {
          ci_data <- tryCatch(
            extract_baf_ci(res, cell, bins),
            error = function(e) { message("  CI extraction failed for ", cell); NULL }
          )
        }
        
        cell_ploidy <- NULL
        if (!is.null(res$allSolutions) && cell %in% names(res$allSolutions)) {
          cell_ploidy <- res$allSolutions[[cell]]$ploidy
        }
        
        profile_plot <- plot_allele_profile(
          nMajor = profiles$nMajor[[cell]],
          nMinor = profiles$nMinor[[cell]],
          bins = bins,
          cell_name = cell,
          ploidy = cell_ploidy,
          segments_ci = ci_data
        )
        
        if (save) {
          ext <- output_format[1]
          profile_file <- file.path(profiles_dir, paste0(cell, ".", ext))
          ggsave(profile_file, profile_plot, width = width, height = height / 2, dpi = dpi, bg = "white")
        }
      }
      message(sprintf("  Saved %d profile plots to: %s", length(cells_to_plot), profiles_dir))
    }

  } else {
    # Non-allele-specific mode (total copy number)
    profiles <- extract_total_profiles(res)
    cell_names <- colnames(profiles$profiles)
    
    # Subset to requested cells if cell_id is provided
    if (!is.null(cell_id)) {
      matching <- cell_id[cell_id %in% cell_names]
      if (length(matching) == 0) {
        warning("No cell_id matched. Available names look like: ",
                paste(head(cell_names, 3), collapse = ", "))
      } else {
        message(sprintf("  Subsetting to %d / %d requested cells",
                        length(matching), length(cell_id)))
        profiles$profiles <- profiles$profiles[, ..matching]
        if (!is.null(profiles$raw)) {
          profiles$raw <- profiles$raw[, ..matching]
        }
        cell_names <- matching
      }
    }
    
    # Create heatmap if requested
    if (plot_type %in% c("all", "heatmap")) {
      message("Creating total copy number heatmap...")
      
      heatmap_plot = plotHeatmap(
        profiles = profiles$profiles,
        bins = bins,
        dendrogram = dendrogram,
        linesize = linesize,
        rownames = show_rownames,
        cell_types = cell_types
      )
      
      if (save) {
        ext <- output_format[1]
        heatmap_file <- file.path(plots_dir, paste0(output_prefix, ".", ext))
        ggsave(heatmap_file, heatmap_plot, width = width, height = height, dpi = dpi, bg = "white")
        message(sprintf("  Saved: %s", heatmap_file))
      }
    }
    
    # Create profile plots if requested
    if (plot_type %in% c("all", "profile")) {
      # Determine which cells to plot
      if (!is.null(cell_id)) {
        cells_to_plot <- cell_id
      } else if (plot_type == "all") {
        cells_to_plot <- cell_names
      } else {
        cells_to_plot <- cell_names[1]
        message(sprintf("No cell_id specified, using first cell: %s", cells_to_plot))
      }
      
      message(sprintf("Creating profile plots for %d cells...", length(cells_to_plot)))
      
      ext <- output_format[1]
      for (cell in cells_to_plot) {
        if (!cell %in% cell_names) {
          warning(sprintf("Cell '%s' not found, skipping", cell))
          next
        }
        
        # Get raw values
        raw_vals <- NULL
        if (!is.null(profiles$raw) && cell %in% colnames(profiles$raw)) {
          raw_vals <- profiles$raw[[cell]]
        }
        
        cell_ploidy <- NULL
        if (!is.null(res$allSolutions) && cell %in% names(res$allSolutions)) {
          cell_ploidy <- res$allSolutions[[cell]]$ploidy
        }
        
        profile_plot = plotProfile(
          segments = profiles$profiles[[cell]],
          raw = raw_vals,
          bins = bins,
          sc = TRUE,
          linesize = linesize,
          cell_name = cell,
          ploidy = cell_ploidy,
          sex = sex
        )
        
        if (save) {
          ext <- output_format[1]
          profile_file <- file.path(profiles_dir, paste0(cell, ".", ext))
          ggsave(profile_file, profile_plot, width = width, height = height / 2, dpi = dpi, bg = "white")
        }
      }
      message(sprintf("  Saved %d profile plots to: %s", length(cells_to_plot), profiles_dir))
    }
  }
  
  message("=== Done ===")
  
  # Return last plot invisibly (or heatmap if "all")
  if (exists("heatmap_plot")) {
    return(invisible(heatmap_plot))
  } else if (exists("profile_plot")) {
    return(invisible(profile_plot))
  }
}

# =============================================================================
# Quick convenience function for common use cases
# =============================================================================

#' Quick plot from RDS file path
#' 
#' @param rds_path Path to ASCAT.sc results RDS file
#' @param ... Additional arguments passed to ascatsc_plot
ascatsc_plot_from_file <- function(rds_path, ...) {
  message(sprintf("Loading: %s", rds_path))
  res <- readRDS(rds_path)
  ascatsc_plot(res, ...)
}
