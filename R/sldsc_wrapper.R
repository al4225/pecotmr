# Stratified LD Score Regression (S-LDSC) post-processing wrappers around polyfun.
#
# This file provides the post-processing layer for the xqtl-protocol sLDSC pipeline:
# read polyfun outputs per trait, compute Gazal-style standardized tau* and
# differential per-SNP heritability (EnrichStat), and run DerSimonian-Laird
# random-effects meta-analysis across traits.
#
# Reference panel convention: all LD-derived quantities (baseline LD scores,
# target LD scores, regression weights, allele frequencies) must come from the
# same reference panel. Do not mix files from different panels (e.g. 1000G vs ADSP).
#
# MAF convention: by default we restrict to MAF > 5% per the sLDSC recommendation.
# Pass maf_cutoff = 0 to opt out (not recommended).
#
# Cross-type comparison: tau* (Gazal et al. 2017 standardization) is the
# cross-type comparable statistic. Use tau* to rank or meta-analyze annotations
# that mix binary and continuous types. E (proportion-based enrichment) is
# scale-dependent for continuous annotations and is only comparable within type.

#' @title Read S-LDSC outputs from polyfun for one trait/run
#'
#' @description Reads the regression outputs produced by `polyfun/ldsc.py` for a
#'   single polyfun run (one trait, one annotation set) and returns them as a
#'   tidy list ready for downstream standardization. Hides the underlying file
#'   formats; downstream code consumes only modeling quantities.
#'
#' @param prefix Character. Path prefix to the polyfun outputs for one trait/run.
#'   The function appends `.results`, `.log`, and `.part_delete` to this prefix.
#'   Example: `"/path/to/cwd/CAD_META.filtered.sumstats.gz"`.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{categories}{Character vector of annotation names, in regression order
#'       (target annotations first, baseline last by convention).}
#'     \item{tau}{Numeric named vector of regression coefficients \eqn{\tau_C}.}
#'     \item{tau_se}{Numeric named vector of standard errors of \eqn{\tau_C}.}
#'     \item{enrichment}{Numeric named vector of \eqn{E_C = \pi^{h^2}_C / \pi^M_C}.}
#'     \item{enrichment_se}{Numeric named vector of standard errors of \eqn{E_C}.}
#'     \item{enrichment_p}{Numeric named vector of p-values for the differential
#'       per-SNP heritability test (the EnrichStat p-value reported by polyfun).}
#'     \item{prop_h2}{Numeric named vector \eqn{\pi^{h^2}_C}.}
#'     \item{prop_snps}{Numeric named vector \eqn{\pi^M_C}.}
#'     \item{h2g}{Scalar total trait heritability \eqn{h^2_g}.}
#'     \item{tau_blocks}{Numeric matrix of per-block jackknife \eqn{\tau} values,
#'       dimensions (n_blocks x n_annotations), columns named by category.}
#'     \item{n_blocks}{Integer number of jackknife blocks (typically 200).}
#'   }
#'
#' @examples
#' \dontrun{
#' run <- read_sldsc_trait("/output/CAD_META.filtered.sumstats.gz")
#' run$tau["my_target_annotation"]
#' }
#'
#' @export
read_sldsc_trait <- function(prefix) {
  stop("Not yet implemented (skeleton).")
}


#' @title Compute per-annotation standard deviation, MAF-restricted
#'
#' @description Computes the standard deviation of each annotation column in the
#'   target annotation files, restricted to SNPs above a MAF cutoff. The MAF
#'   restriction is required for the standardization to be internally consistent
#'   with polyfun's regression, which operates on MAF > cutoff SNPs by default.
#'
#' @param target_anno_dir Character. Directory containing target annotation files
#'   (one per chromosome) in polyfun's `.annot.gz` format.
#' @param frqfile_dir Character or NULL. Directory containing PLINK `.frq` files
#'   for the reference panel. Required when `maf_cutoff > 0`; the function
#'   errors if missing.
#' @param plink_name Character. Filename prefix of the `.frq` files
#'   (e.g. `"ADSP_chr"`). Files are expected at `{plink_name}{chr}.frq`.
#' @param maf_cutoff Numeric, default `0.05`. Only SNPs with `MAF > maf_cutoff`
#'   contribute to sd. Set to `0` to disable filtering (must match the
#'   `--not-M-5-50` regression mode for internal consistency).
#' @param annot_cols Character or integer vector, default NULL. Annotation columns
#'   to compute sd for. If NULL, all annotation columns (i.e. all columns past
#'   the standard fixed columns CHR/SNP/BP/CM/A1/A2/MAF) are used.
#'
#' @return Named numeric vector of \eqn{sd_C} values, one per annotation.
#'
#' @examples
#' \dontrun{
#' sd_annot <- compute_sldsc_annot_sd(
#'   target_anno_dir = "/output/ldscores/AC_DeJager_eQTL",
#'   frqfile_dir     = "/ref/ADSP/frq",
#'   plink_name      = "ADSP_chr",
#'   maf_cutoff      = 0.05
#' )
#' }
#'
#' @export
compute_sldsc_annot_sd <- function(target_anno_dir, frqfile_dir = NULL,
                                   plink_name = "ADSP_chr",
                                   maf_cutoff = 0.05, annot_cols = NULL) {
  stop("Not yet implemented (skeleton).")
}


#' @title Reference-panel SNP count at a given MAF cutoff
#'
#' @description Returns the number of SNPs in the reference panel above the MAF
#'   cutoff, matching the regression's M_5_50 (when `maf_cutoff > 0`) or M
#'   (when `maf_cutoff == 0`). Used as \eqn{M_{ref}} in the standardization.
#'
#' @param target_anno_dir Character. Directory containing per-chromosome
#'   `.l2.M_5_50` (when MAF-filtered) or `.l2.M` (when not) files produced by
#'   polyfun's `compute_ldscores.py` for the target annotation.
#' @param maf_cutoff Numeric, default `0.05`.
#'
#' @return Scalar integer: total number of SNPs in the reference panel above the cutoff.
#'
#' @export
compute_sldsc_M_ref <- function(target_anno_dir, maf_cutoff = 0.05) {
  stop("Not yet implemented (skeleton).")
}


#' @title Detect whether each annotation is binary or continuous
#'
#' @description Inspects each annotation column in the target annotation files
#'   and returns whether its values lie in \{0, 1\} (binary) or take other values
#'   (continuous). Used to select the appropriate within-type headline statistic.
#'   For cross-type comparison, always use \eqn{\tau^*_C} regardless of the flag.
#'
#' @param target_anno_dir Character. Directory containing the target `.annot.gz`
#'   files (one per chromosome).
#' @param annot_cols Character or integer vector, default NULL.
#'
#' @return Named logical vector: TRUE for binary, FALSE for continuous.
#'
#' @export
is_binary_sldsc_annot <- function(target_anno_dir, annot_cols = NULL) {
  stop("Not yet implemented (skeleton).")
}


#' @title Standardize tau and compute per-block tau* and EnrichStat for one polyfun run
#'
#' @description Given the polyfun read result plus the annotation sd's and
#'   reference SNP count, computes the Gazal-style standardized tau (\eqn{\tau^*})
#'   and the differential per-SNP heritability statistic (EnrichStat) for the
#'   target annotations, including their per-block jackknife values for use in
#'   cross-trait meta-analysis.
#'
#' @details The standardization is
#'   \deqn{\tau^*_C = \tau_C \cdot sd_C \cdot M_{ref} / h^2_g}
#'   applied both to the point estimate and to each of the n_blocks jackknife
#'   blocks. The EnrichStat is
#'   \deqn{\frac{h^2_g}{M_{ref}}\left[\frac{\pi^{h^2}_C}{\pi^M_C}
#'         - \frac{1-\pi^{h^2}_C}{1-\pi^M_C}\right]}
#'   computed from the per-block tau values, with jackknife SE from the per-block
#'   variance \eqn{\sqrt{\frac{(B-1)^2}{B} Var_b}}.
#'
#'   Only target annotations are returned; baseline annotations are dropped after
#'   serving their role of conditioning the regression.
#'
#' @param trait_data List. Output of \code{\link{read_sldsc_trait}} for one polyfun run.
#' @param sd_annot Named numeric vector. Output of
#'   \code{\link{compute_sldsc_annot_sd}}, restricted to target annotations.
#' @param M_ref Scalar. Output of \code{\link{compute_sldsc_M_ref}}.
#' @param target_categories Character vector or NULL. Annotation names to keep.
#'   If NULL, all categories with names matching `names(sd_annot)` are kept.
#'
#' @return A list with components:
#'   \describe{
#'     \item{summary}{Data frame with one row per target annotation and columns:
#'       `target`, `tau`, `tau_se`, `tau_star`, `tau_star_se`, `enrichment`,
#'       `enrichment_se`, `enrichment_p`, `enrichstat`, `enrichstat_se`.}
#'     \item{tau_star_blocks}{Matrix (n_blocks x n_target) of per-block \eqn{\tau^*_C}.}
#'     \item{enrichstat_blocks}{Matrix (n_blocks x n_target) of per-block EnrichStat.}
#'     \item{h2g}{Scalar trait heritability.}
#'     \item{n_blocks}{Integer.}
#'   }
#'
#' @export
standardize_sldsc_trait <- function(trait_data, sd_annot, M_ref,
                                    target_categories = NULL) {
  stop("Not yet implemented (skeleton).")
}


#' @title Random-effects meta-analysis of S-LDSC quantities across traits
#'
#' @description DerSimonian-Laird random-effects meta-analysis of one S-LDSC
#'   quantity across multiple traits, for one annotation. Used as the entry
#'   point both for the default pipeline meta and for user-driven re-meta over
#'   subsets of traits.
#'
#' @details Implements
#'   \deqn{\hat\theta_{meta} = \sum_i w_i \hat\theta_i / \sum_i w_i,
#'         \quad SE_{meta} = 1/\sqrt{\sum_i w_i},
#'         \quad w_i = 1/(SE_i^2 + \hat\sigma^2)}
#'   via `rmeta::meta.summaries(..., method = "random")`. The two-sided p-value
#'   is computed from the meta z-score under the standard normal.
#'
#'   `quantity = "tau_star"` and `"enrichstat"` use the jackknife SE from
#'   per-block delete values (computed inside
#'   \code{\link{standardize_sldsc_trait}}). `quantity = "enrichment"` uses
#'   the SE reported directly by the regression engine.
#'
#' @param per_trait_estimates Named list of standardized per-trait results from
#'   \code{\link{standardize_sldsc_trait}}. Pass a subset to re-meta on a
#'   user-chosen group of traits.
#' @param category Character. Annotation name to meta-analyze.
#' @param quantity Character, one of `"tau_star"`, `"enrichment"`, or
#'   `"enrichstat"`.
#'
#' @return A list with components: `mean`, `se`, `p`, `n_traits`,
#'   `traits_used`, and `tau2` (the DerSimonian-Laird between-trait variance).
#'
#' @examples
#' \dontrun{
#' meta_all <- meta_sldsc_random(per_trait_results, "my_anno", "tau_star")
#'
#' # subset to a custom trait group
#' subset <- per_trait_results[c("CAD_META", "AD_GWAX", "PD_meta")]
#' meta_neuro <- meta_sldsc_random(subset, "my_anno", "tau_star")
#' }
#'
#' @export
meta_sldsc_random <- function(per_trait_estimates, category,
                              quantity = c("tau_star", "enrichment", "enrichstat")) {
  stop("Not yet implemented (skeleton).")
}


#' @title End-to-end S-LDSC post-processing across traits, single + joint in one pass
#'
#' @description Top-level convenience wrapper. Reads polyfun outputs (single-tau
#'   runs and the joint-tau run) for each trait, standardizes both modes, and
#'   runs the default random-effects meta-analysis across all traits supplied.
#'   Single and joint analyses are produced in one pass; no `joint_tau` flag.
#'
#' @details For a set of \eqn{N} target annotations and one trait, the caller is
#'   expected to have already produced \eqn{N+1} polyfun runs: \eqn{N} single-tau
#'   runs (each fitting one target + baseline) and 1 joint-tau run (fitting all
#'   \eqn{N} targets + baseline together). The pipeline reads all of these, drops
#'   baseline categories, standardizes both modes, and assembles consolidated
#'   per-trait and meta data frames.
#'
#'   Reports baseline annotation count and names via `message()` once at the
#'   start of the run for transparency.
#'
#'   Cross-type comparison: the `meta$tau_star` frame is the apple-to-apple table
#'   for ranking annotations regardless of binary/continuous type. The `is_binary`
#'   column lets callers filter further when within-type comparison is desired.
#'
#' @param trait_single_prefixes Named list. For each trait (name = trait id),
#'   a character vector of length \eqn{N} giving the polyfun output prefixes for
#'   the \eqn{N} single-tau runs (one per target annotation). Order must match
#'   `target_categories`.
#' @param trait_joint_prefix Named character. For each trait (name = trait id,
#'   matching `trait_single_prefixes`), the polyfun output prefix for the joint-tau
#'   run that fits all \eqn{N} targets together.
#' @param target_anno_dir Character. Directory of the target annotation files
#'   (used for `sd_C` and binary detection). When all \eqn{N} targets share
#'   one directory (multi-column annot files), pass that. When each target has
#'   its own directory, pass a named character vector (names = target categories).
#' @param frqfile_dir Character or NULL. Directory of `.frq` files for the panel.
#' @param plink_name Character. Filename prefix of `.frq` files (default `"ADSP_chr"`).
#' @param maf_cutoff Numeric, default `0.05`.
#' @param target_categories Character vector or NULL. Target annotation names to
#'   report. If NULL, auto-detected from the joint-run results as all categories
#'   not in baseline.
#'
#' @return A list with components:
#'   \describe{
#'     \item{per_trait}{Named list of per-trait results. For each trait, a list with:
#'       `summary` (wide data frame of single + joint estimates side-by-side per
#'       target annotation, with `is_binary` flag), `tau_star_blocks_single`,
#'       `tau_star_blocks_joint`, `enrichstat_blocks_single`,
#'       `enrichstat_blocks_joint` (all matrices), and `h2g`.}
#'     \item{meta}{List of three data frames, each with one row per target annotation:
#'       \itemize{
#'         \item `tau_star`: columns `target`, `is_binary`, `single_mean`,
#'           `single_se`, `single_p`, `joint_mean`, `joint_se`, `joint_p`, `n_traits`.
#'         \item `enrichment`: same shape.
#'         \item `enrichstat`: same shape.
#'       }
#'       The `tau_star` frame is the cross-type comparable headline.}
#'     \item{params}{List echoing `maf_cutoff`, `M_ref`, `target_categories`,
#'       `n_baseline`, `baseline_categories`, and `trait_names` for reproducibility.}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- sldsc_postprocessing_pipeline(
#'   trait_single_prefixes = list(
#'     CAD_META = c("/output/CAD_META.single_anno1", "/output/CAD_META.single_anno2"),
#'     AD_GWAX  = c("/output/AD_GWAX.single_anno1",  "/output/AD_GWAX.single_anno2")
#'   ),
#'   trait_joint_prefix = c(
#'     CAD_META = "/output/CAD_META.joint",
#'     AD_GWAX  = "/output/AD_GWAX.joint"
#'   ),
#'   target_anno_dir = "/output/ldscores/my_targets",
#'   frqfile_dir     = "/ref/ADSP/frq",
#'   plink_name      = "ADSP_chr",
#'   maf_cutoff      = 0.05,
#'   target_categories = c("anno1", "anno2")
#' )
#' res$meta$tau_star      # cross-type comparable headline
#' res$meta$enrichment    # within-binary headline
#'
#' # later, re-meta on a custom subset of traits
#' meta_neuro <- meta_sldsc_random(
#'   res$per_trait[c("AD_GWAX", "PD_meta")], "anno1", "tau_star"
#' )
#' }
#'
#' @seealso \code{\link{read_sldsc_trait}}, \code{\link{standardize_sldsc_trait}},
#'   \code{\link{meta_sldsc_random}}
#'
#' @export
sldsc_postprocessing_pipeline <- function(trait_single_prefixes,
                                          trait_joint_prefix,
                                          target_anno_dir,
                                          frqfile_dir = NULL,
                                          plink_name = "ADSP_chr",
                                          maf_cutoff = 0.05,
                                          target_categories = NULL) {
  stop("Not yet implemented (skeleton).")
}
