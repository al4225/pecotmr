#' Load a PLINK .bim file for cTWAS
#'
#' @description
#' \strong{Deprecated.} Use [read_bim()] via the standard I/O path
#' instead. This wrapper remains for backwards compatibility and calls
#' [read_bim()] internally, mapping its output to the legacy column names.
#'
#' @param bim_file_path Path to a PLINK \code{.bim} file (or a \code{.bed}
#'   file - the \code{.bim} extension is resolved automatically).
#'
#' @return A data.frame with columns \code{chrom}, \code{id}, \code{GD},
#'   \code{pos}, \code{A1}, \code{A2}. Variant IDs are normalised via
#'   [normalize_variant_id()].
#'
#' @export
ctwas_bimfile_loader <- function(bim_file_path) {
  .Deprecated("read_bim", package = "pecotmr",
              msg = "ctwas_bimfile_loader() is deprecated. Use read_bim() instead.")
  # read_bim() expects a .bed path and derives .bim from it.
  # Accept either .bim or .bed and normalise to .bed.
  bed_path <- sub("\\.bim$", ".bed", bim_file_path)
  bim <- read_bim(bed_path)
  # Map new column names back to legacy names
  snp_info <- data.frame(
    chrom = bim$chrom,
    id    = normalize_variant_id(bim$id),
    GD    = bim$gpos,
    pos   = bim$pos,
    A1    = bim$a1,
    A2    = bim$a0,
    stringsAsFactors = FALSE
  )
  return(snp_info)
}

#' Load cTWAS LD meta-data
#'
#' @description
#' \strong{Deprecated.} Use [ld_loader()] with its \code{LD_info}
#' argument instead. This wrapper remains for backwards compatibility and
#' produces the same \code{list(LD_info, region_info)} output as the original.
#'
#' @param ld_meta_data_file Path to the LD meta-data TSV file.
#' @param subset_region_ids Optional character vector of region IDs
#'   (\code{"chrom_start_end"}) to subset to.
#'
#' @return A list with components:
#' \describe{
#'   \item{LD_info}{Data.frame with columns \code{region_id}, \code{LD_file},
#'     \code{SNP_file}.}
#'   \item{region_info}{Data.frame with columns \code{chrom}, \code{start},
#'     \code{stop}, \code{region_id}.}
#' }
#'
#' @importFrom vroom vroom
#' @export
get_ctwas_meta_data <- function(ld_meta_data_file, subset_region_ids = NULL) {
  .Deprecated("ld_loader", package = "pecotmr",
              msg = "get_ctwas_meta_data() is deprecated. Use ld_loader() with LD_info instead.")
  LD_info <- as.data.frame(vroom(ld_meta_data_file))
  colnames(LD_info)[1] <- "chrom"
  LD_info$region_id <- paste(as.integer(strip_chr_prefix(LD_info$chrom)),
                             LD_info$start, LD_info$end, sep = "_")
  LD_info$LD_file <- paste0(dirname(ld_meta_data_file), "/",
                            gsub(",.*$", "", LD_info$path))
  LD_info$SNP_file <- paste0(LD_info$LD_file, ".bim")
  LD_info <- LD_info[, c("region_id", "LD_file", "SNP_file")]
  region_info <- LD_info[, "region_id", drop = FALSE]
  region_info$chrom <- as.integer(gsub("\\_.*$", "", region_info$region_id))
  region_info$start <- as.integer(gsub("\\_.*$", "",
                                       sub("^.*?\\_", "", region_info$region_id)))
  region_info$stop <- as.integer(sub("^.*?\\_", "",
                                      sub("^.*?\\_", "", region_info$region_id)))
  region_info$region_id <- paste0(region_info$chrom, "_",
                                   region_info$start, "_",
                                   region_info$stop)
  region_info <- region_info[, c("chrom", "start", "stop", "region_id")]
  if (!is.null(subset_region_ids)) {
    region_info <- region_info[region_info$region_id %in% subset_region_ids, ]
  }
  return(list(LD_info = LD_info, region_info = region_info))
}

#' Function to select variants for ctwas weights input
#' @param region_data A list of list containing weights list and snp_info list data for multiple genes/events within a single LD block region.
#' @param export_twas_weight_db A list of list of fine-mapping result data formatted by generate_twas_db function.
#' @param region_block A string for region information for region_weights, consisted of chromosome number, star and end position of LD block conneced with "_".
#' @export
trim_ctwas_variants <- function(region_data, twas_weight_cutoff = 1e-5, cs_min_cor = 0.8,
                                min_pip_cutoff = 0.1, max_num_variants = 1000) {
  # internal functions to select variants for a gene-context pair weight list
  select_variants <- function(group_name, study_name, region_data, cs_min_cor, min_pip_cutoff, max_num_variants) {
    weight_list <- region_data$weights[[group_name]][[study_name]]
    context <- weight_list$context
    selected_variants_by_context <- c()
    molecular_id <- gsub("\\|.*", "", group_name)

    if ("cs_variants" %in% names(region_data$susie_weights_intermediate[[molecular_id]][[context]]) & length(region_data$susie_weights_intermediate[[molecular_id]][[context]][["cs_variants"]]) != 0) {
      cs_min_abs_cor <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_purity$min.abs.corr
      for (L in seq_along(region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants)) {
        # we includ all variants in $cs_variant if min_abs_corr > cs_min_cor for the set
        if (cs_min_abs_cor[L] >= cs_min_cor) {
          cs_variants <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants[[L]]
          selected_variants_by_context <- cs_variants[cs_variants %in% rownames(weight_list$wgt)]
        }
      }
    }
    context_pip <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$pip
    # variant IDs are in canonical chr-prefix format from allele_qc pipeline
    high_pip_variants <- names(context_pip[context_pip > min_pip_cutoff])[names(context_pip[context_pip > min_pip_cutoff]) %in% rownames(weight_list$wgt)]
    selected_variants_by_context <- unique(c(selected_variants_by_context, high_pip_variants))

    # prioritize SNPs based on PIP if max_num_variants different from Inf
    available_variants <- intersect(rownames(weight_list$wgt), names(context_pip))
    prioritized <- unique(c(selected_variants_by_context, setdiff(available_variants, selected_variants_by_context)))
    prioritized <- prioritized[order(-context_pip[prioritized])]
    selected_variants_by_context <- head(prioritized, max_num_variants)
    weight_list$wgt <- weight_list$wgt[selected_variants_by_context, , drop = FALSE]
    return(weight_list)
  }
  merge_by_study <- function(weights) {
    weight_list <- list()
    for (group in names(weights)) {
      for (study in names(weights[[group]])) {
        weight_list[[study]][[group]] <- weights[[group]][[study]]
      }
    }
    return(weight_list)
  }

  weights <- setNames(lapply(names(region_data$weights), function(group) {
    for (study in names(region_data$weights[[group]])) {
      region_data$weights[[group]][[study]]$wgt <- region_data$weights[[group]][[study]]$wgt[abs(region_data$weights[[group]][[study]]$wgt[, 1]) >= twas_weight_cutoff, , drop = FALSE]
      if (nrow(region_data$weights[[group]][[study]]$wgt) < 1) {
        region_data$weights[[group]][[study]] <- NULL
        next
      }
      if (all(is.na(region_data$weights[[group]][[study]]$wgt[, 1])) || all(is.nan(region_data$weights[[group]][[study]]$wgt[, 1]))) {
        region_data$weights[[group]][[study]] <- NULL
        next
      }
      if (nrow(region_data$weights[[group]][[study]]$wgt) < max_num_variants) {
        region_data$weights[[group]][[study]]$n_wgt <- nrow(region_data$weights[[group]][[study]]$wgt)
      } else {
        region_data$weights[[group]][[study]] <- select_variants(group, study, region_data, cs_min_cor = cs_min_cor, min_pip_cutoff = min_pip_cutoff, max_num_variants = max_num_variants)
        region_data$weights[[group]][[study]]$n_wgt <- nrow(region_data$weights[[group]][[study]]$wgt)
      }
      region_data$weights[[group]] <- Filter(Negate(is.null), region_data$weights[[group]])
      context_range <- as.integer(sapply(rownames(region_data$weights[[group]][[study]]$wgt), function(variant) strsplit(variant, "\\:")[[1]][2]))
      if(twas_weight_cutoff!=0 | cs_min_cor!=0 | min_pip_cutoff!=0 | max_num_variants!=Inf){
        region_data$weights[[group]][[study]][["p0"]] = min(context_range)# update min max position
        region_data$weights[[group]][[study]][["p1"]] = max(context_range)
      }
    }
    return(region_data$weights[[group]])
  }), names(region_data$weights))
  weights <- Filter(Negate(is.null), weights)
  weights <- merge_by_study(weights)
  return(weights)
}
