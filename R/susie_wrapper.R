#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param prior_weights Optional numeric vector of prior weights for each element in lbf.
#' @return A named numeric vector of alpha values corresponding to the input LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbf_to_alpha_vector(lbf)
#' print(alpha)
#' @noRd
lbf_to_alpha_vector <- function(lbf, prior_weights = NULL) {
  if (is.null(prior_weights)) prior_weights <- rep(1 / length(lbf), length(lbf))
  maxlbf <- max(lbf)

  # If maxlbf is 0, return a vector of zeros
  if (maxlbf == 0) {
    return(setNames(rep(0, length(lbf)), names(lbf)))
  }

  # w is proportional to BF, subtract max for numerical stability
  w <- exp(lbf - maxlbf)

  # Posterior prob for each SNP
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  return(alpha)
}

#' Applies the 'lbf_to_alpha_vector' function row-wise to a matrix of log Bayes factors
#' to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF matrix.
#' @examples
#' lbf_matrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alpha_matrix <- lbf_to_alpha(lbf_matrix)
#' print(alpha_matrix)
#' @export
lbf_to_alpha <- function(lbf) {
  alpha_matrix <- t(apply(as.matrix(lbf), 1, lbf_to_alpha_vector))
  if (ncol(lbf) == 1) alpha_matrix <- matrix(alpha_matrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
  return(alpha_matrix)
}

format_pip_column <- function(method) {
  paste0("pip_", method)
}

resolve_pip_column <- function(top_loci, method = NULL) {
  if (is.null(top_loci) || nrow(top_loci) == 0) return(NULL)
  if (!is.null(method)) {
    pip_col <- format_pip_column(method)
    if (pip_col %in% names(top_loci)) return(pip_col)
  }
  if ("pip" %in% names(top_loci)) return("pip")
  pip_cols <- grep("^pip_", names(top_loci), value = TRUE)
  if (length(pip_cols) == 1) return(pip_cols)
  NULL
}

format_cs_column <- function(coverage, method) {
  pct <- as.numeric(coverage) * 100
  if (is.na(pct)) stop("coverage must be numeric.")
  label <- if (abs(pct - round(pct)) < 1e-8) {
    as.character(as.integer(round(pct)))
  } else {
    gsub("\\.", "_", format(pct, scientific = FALSE, trim = TRUE))
  }
  paste0("CS_", label, "_", method)
}

.translate_legacy_cs_column_name <- function(coverage) {
  if (is.null(coverage)) return(NULL)
  vapply(coverage, function(x) {
    x <- as.character(x)
    old_match <- regexec("^cs_coverage_([0-9.]+)$", x, ignore.case = TRUE)
    old_parts <- regmatches(x, old_match)[[1]]
    if (length(old_parts) == 2) return(format_cs_column(as.numeric(old_parts[[2]]), "susie"))
    x
  }, character(1), USE.NAMES = FALSE)
}

.translate_legacy_top_loci_cs_columns <- function(top_loci) {
  if (!is.data.frame(top_loci)) return(top_loci)
  names(top_loci) <- .translate_legacy_cs_column_name(names(top_loci))
  if ("pip_susie" %in% names(top_loci) && !"pip" %in% names(top_loci)) {
    names(top_loci)[names(top_loci) == "pip_susie"] <- "pip"
  }
  top_loci
}

.set_finemapping_fit_class <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  method_class <- switch(method,
    susie = "susie",
    susie_inf = "susie_inf",
    susie_rss = "susie_rss",
    single_effect = "susie_rss",
    bayesian_conditional_regression = "susie_rss",
    fsusie = "susiF",
    mvsusie = "mvsusie",
    NULL
  )
  if (!is.null(method_class)) class(fit) <- unique(c(method_class, class(fit)))
  fit
}

prepare_susie_from_inf_args <- function(args, susie_inf_fit, refine_default = NULL) {
  L <- args[["L"]]
  if (is.null(L)) L <- length(susie_inf_fit$V)
  if (is.null(args[["refine"]]) && !is.null(refine_default)) args[["refine"]] <- refine_default
  args[["unmappable_effects"]] <- "none"
  args[["model_init"]] <- susie_inf_fit
  if (!is.null(args[["L_greedy"]])) args[["L_greedy"]] <- min(length(susie_inf_fit$V), L)
  args
}

fit_susie_inf_then_susie <- function(X, y, args = list(),
                                     susie_inf_args = list(),
                                     susie_args = list(),
                                     fitted_models = NULL) {
  if (is.null(fitted_models)) fitted_models <- list()
  susie_inf_fit <- fitted_models[["susie_inf"]]
  susie_fit <- fitted_models[["susie"]]

  if (is.null(susie_inf_fit)) {
    fit_args <- modifyList(args, susie_inf_args)
    fit_args <- modifyList(fit_args, list(
      X = X, y = y, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susie_inf_fit <- do.call(susie, fit_args)
  }
  susie_inf_fit <- .set_finemapping_fit_class(susie_inf_fit, "susie_inf")

  if (is.null(susie_fit)) {
    fit_args <- prepare_susie_from_inf_args(modifyList(args, susie_args), susie_inf_fit, refine_default = TRUE)
    susie_fit <- do.call(susie, c(list(X = X, y = y), fit_args))
  }
  susie_fit <- .set_finemapping_fit_class(susie_fit, "susie")

  list(susie = susie_fit, susie_inf = susie_inf_fit)
}

#' Post-process Fine-mapping Fits
#'
#' Applies method-aware post-processing to one or more SuSiE-family fits and
#' builds both a method-specific result list and shared top-loci tables.
#'
#' @param fits Named list of fine-mapping fits. Names define method identity,
#'   for example \code{susie}, \code{susie_inf}, \code{susie_rss},
#'   \code{mvsusie}, or \code{fsusie}.
#' @param data_x Genotype matrix, LD/correlation matrix, or other method-specific
#'   input used for credible-set purity and correlations.
#' @param data_y Phenotype vector/matrix or summary statistics. Default NULL.
#' @param X_scalar Scaling factor for genotype effects. Default 1.
#' @param y_scalar Scaling factor for phenotype effects. Default 1.
#' @param maf Minor allele frequencies. Default NULL.
#' @param coverage Primary credible-set coverage.
#' @param secondary_coverage Additional credible-set coverages.
#' @param signal_cutoff PIP cutoff for including non-CS variants in top loci.
#' @param other_quantities Optional list carried into each method result.
#' @param prior_eff_tol Tolerance for retaining effects by prior variance.
#' @param min_abs_corr Minimum absolute correlation for credible-set purity.
#' @return A list with \code{finemapping_results} (per-method post-processed
#'   objects, each carrying a trimmed fit and method-specific intermediates)
#'   and the single unified \code{top_loci} table in the fixed 22-column
#'   shape (see \code{\link{build_top_loci}}). Per-method contributions are
#'   produced by \code{build_top_loci()} once per method and row-bound into
#'   this single \code{top_loci}. There is no separately exposed
#'   \code{top_loci_long} or wide-format \code{top_loci}.
#' @export
postprocess_finemapping_fits <- function(fits, data_x, data_y = NULL,
                                         X_scalar = 1, y_scalar = 1,
                                         maf = NULL, coverage = NULL,
                                         secondary_coverage = c(0.7, 0.5),
                                         signal_cutoff = 0.1,
                                         other_quantities = NULL,
                                         prior_eff_tol = 1e-9,
                                         min_abs_corr = 0.8) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) == 0) stop("At least one fine-mapping fit must be supplied.")
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("fits must be a named list; names define method identity.")
  }

  # One method for-loop: each method calls build_top_loci() once per fit; the
  # per-method 22-column contributions are row-bound below into the single
  # final `top_loci` table. There is no separately exposed long or wide table.
  posts <- lapply(names(fits), function(method) {
    fit <- .set_finemapping_fit_class(fits[[method]], method)
    postprocess_finemapping_fit(
      fit, method = method, data_x = data_x, data_y = data_y,
      X_scalar = X_scalar, y_scalar = y_scalar, maf = maf,
      coverage = coverage, secondary_coverage = secondary_coverage,
      signal_cutoff = signal_cutoff, other_quantities = other_quantities,
      prior_eff_tol = prior_eff_tol, min_abs_corr = min_abs_corr
    )
  })
  names(posts) <- names(fits)

  per_method <- lapply(posts, function(x) x$top_loci)
  per_method <- per_method[!vapply(per_method, is.null, logical(1))]
  top_loci <- if (length(per_method) == 0L) {
    .empty_top_loci()
  } else {
    do.call(rbind, per_method)
  }
  rownames(top_loci) <- NULL
  posts <- lapply(posts, function(x) {
    x$top_loci <- NULL
    x
  })

  list(
    finemapping_results = posts,
    top_loci = top_loci
  )
}

postprocess_finemapping_fit <- function(fit, ...) {
  UseMethod("postprocess_finemapping_fit")
}

#' @exportS3Method
postprocess_finemapping_fit.susie <- function(fit, method = "susie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_inf <- function(fit, method = "susie_inf", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susie_rss <- function(fit, method = "susie_rss", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "Xcorr", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.mvsusie <- function(fit, method = "mvsusie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "X", ...)
}

#' @exportS3Method
postprocess_finemapping_fit.susiF <- function(fit, method = "fsusie", ...) {
  .postprocess_finemapping_fit_common(fit, method = method, cs_input = "fsusie", ...)
}

.postprocess_finemapping_fit_common <- function(fit, method, data_x, data_y = NULL,
                                                X_scalar = 1, y_scalar = 1,
                                                maf = NULL, coverage = NULL,
                                                secondary_coverage = c(0.7, 0.5),
                                                signal_cutoff = 0.1,
                                                other_quantities = NULL,
                                                prior_eff_tol = 1e-9,
                                                min_abs_corr = 0.8,
                                                cs_input = c("X", "Xcorr", "fsusie")) {
  cs_input <- match.arg(cs_input)
  variant_names <- extract_variant_names(fit)
  sumstats <- extract_sumstats(fit, data_x, data_y, X_scalar, y_scalar, method)
  effect_idx <- select_effects(fit, prior_eff_tol)
  cs_tables <- compute_cs_tables(
    fit, data_x = data_x, coverage = coverage,
    secondary_coverage = secondary_coverage, method = method,
    cs_input = cs_input, min_abs_corr = min_abs_corr
  )
  top_loci <- build_top_loci(
    fit, cs_tables, variant_names = variant_names, sumstats = sumstats,
    maf = maf, method = method, signal_cutoff = signal_cutoff,
    data_x = data_x, data_y = data_y, other_quantities = other_quantities
  )

  trimmed <- trim_finemapping_fit(fit, effect_idx, method, cs_tables)

  # Build FineMappingResult S4 object. The S4 contract (validity check,
  # vcf_writer, getPIP, getCS) still expects `variant_id`, `pip`, and an
  # integer `cs` column on the slot. To avoid rippling renames into
  # AllClasses / AllMethods / vcf_writer for this change, we project the
  # new 22-column `top_loci` into the legacy slot shape here, in
  # susie_wrapper only. The wrapper-facing `top_loci` returned to callers
  # is unchanged.
  s4_top_loci <- .top_loci_for_s4_slot(top_loci)
  fm_result <- FineMappingResult(
    variant_names = variant_names,
    trimmed_fit = trimmed,
    top_loci = s4_top_loci,
    method = method,
    sumstats = sumstats
  )

  # Also return as list for backwards compatibility with existing consumers
  res <- list(
    variant_names = variant_names,
    result_trimmed = trimmed,
    top_loci = top_loci,
    finemapping_result = fm_result
  )
  if (!is.null(sumstats)) res$sumstats <- sumstats
  sample_names <- .sample_names_from_data_y(data_y)
  if (!is.null(sample_names)) res$sample_names <- sample_names
  if (method == "mvsusie" && !is.null(fit$outcome_names)) res$context_names <- fit$outcome_names
  analysis_script <- load_script()
  if (analysis_script != "") res$analysis_script <- analysis_script
  if (!is.null(other_quantities)) res$other_quantities <- other_quantities
  res
}

extract_variant_names <- function(fit) {
  variant_names <- names(fit$pip)
  if (is.null(variant_names)) variant_names <- colnames(fit$alpha)
  if (is.null(variant_names)) variant_names <- paste0("variant_", seq_along(fit$pip))
  tryCatch(normalize_variant_id(variant_names), error = function(e) variant_names)
}

extract_sumstats <- function(fit, data_x, data_y, X_scalar = 1, y_scalar = 1, method = "susie") {
  if (is.null(data_y)) return(NULL)
  if (method == "susie_rss") return(data_y)
  if (is.list(data_y) && !is.data.frame(data_y) &&
      any(c("betahat", "sebetahat", "z") %in% names(data_y))) {
    return(data_y)
  }
  if (is.null(data_x)) return(NULL)
  if (is.matrix(data_y) || is.data.frame(data_y)) {
    if (ncol(as.matrix(data_y)) != 1) return(NULL)
  }
  sumstats <- univariate_regression(data_x, data_y)
  y_scalar <- if (is.null(y_scalar) || all(y_scalar == 1)) 1 else y_scalar
  X_scalar <- if (is.null(X_scalar) || all(X_scalar == 1)) 1 else X_scalar
  sumstats$betahat <- sumstats$betahat * y_scalar / X_scalar
  sumstats$sebetahat <- sumstats$sebetahat * y_scalar / X_scalar
  sumstats
}

.sample_names_from_data_y <- function(data_y) {
  if (is.null(data_y) || is.list(data_y)) return(NULL)
  rownames(as.matrix(data_y))
}

select_effects <- function(fit, prior_eff_tol = 1e-9) {
  alpha <- .as_effect_matrix(fit$alpha)
  n_effects <- nrow(alpha)
  if (n_effects == 0) return(integer(0))
  if (!is.null(fit$V)) {
    which(fit$V > prior_eff_tol)
  } else {
    seq_len(n_effects)
  }
}

.as_effect_matrix <- function(x) {
  if (is.null(x)) return(matrix(numeric(0), nrow = 0))
  if (is.list(x) && !is.data.frame(x)) return(do.call(rbind, x))
  as.matrix(x)
}

.as_lbf_matrix <- function(fit) {
  if (!is.null(fit$lbf_variable)) return(.as_effect_matrix(fit$lbf_variable))
  if (!is.null(fit$lBF)) return(.as_effect_matrix(fit$lBF))
  NULL
}

#' @importFrom susieR get_cs_correlation
#' @noRd
compute_cs_tables <- function(fit, data_x, coverage = NULL,
                              secondary_coverage = c(0.7, 0.5),
                              method = "susie", cs_input = c("X", "Xcorr", "fsusie"),
                              min_abs_corr = 0.8) {
  cs_input <- match.arg(cs_input)
  primary_coverage <- coverage
  if (is.null(primary_coverage)) primary_coverage <- fit$sets$requested_coverage
  if (is.null(primary_coverage)) primary_coverage <- 0.95
  coverages <- unique(c(primary_coverage, secondary_coverage))
  coverages <- coverages[!is.na(coverages)]

  tables <- lapply(coverages, function(cov) {
    compute_cs_table(fit, data_x, coverage = cov, cs_input = cs_input, min_abs_corr = min_abs_corr)
  })
  names(tables) <- vapply(coverages, format_cs_column, character(1), method = method)
  attr(tables, "coverage") <- coverages
  tables
}

compute_cs_table <- function(fit, data_x, coverage, cs_input = c("X", "Xcorr", "fsusie"),
                             min_abs_corr = 0.8) {
  cs_input <- match.arg(cs_input)
  if (cs_input == "fsusie") {
    sets <- tryCatch(
      fsusie_get_cs(fit, data_x, requested_coverage = coverage),
      error = function(e) list(cs = list(), requested_coverage = coverage)
    )
    if (is.null(sets$cs) || length(sets$cs) == 0 || all(vapply(sets$cs, is.null, logical(1)))) {
      sets$cs <- list()
      return(list(sets = sets, cs_corr = NULL, pip = fit$pip))
    }
    tmp <- fit
    tmp$sets <- sets
    cs_corr <- if (requireNamespace("fsusieR", quietly = TRUE)) {
      tryCatch(fsusieR::cal_cor_cs(tmp, data_x), error = function(e) NULL)
    } else {
      NULL
    }
    return(list(sets = sets, cs_corr = cs_corr, pip = fit$pip))
  }

  if (cs_input == "X") {
    sets <- susie_get_cs(fit, X = data_x, coverage = coverage, min_abs_corr = min_abs_corr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, X = data_x)
  } else {
    sets <- susie_get_cs(fit, Xcorr = data_x, coverage = coverage, min_abs_corr = min_abs_corr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, Xcorr = data_x)
  }
  out
}

#' Build the unified single top-loci table for one fit and one method.
#'
#' Produces the per-fit, per-method contribution to the unified \code{top_loci}
#' table in the fixed 22-column shape. The outer
#' \code{postprocess_finemapping_fits()} loop calls this once per method per fit
#' and row-binds the results into the single final \code{top_loci} table that is
#' exposed by \code{format_finemapping_output()}.
#'
#' This function replaces the previous \code{build_top_loci_long()},
#' \code{build_top_loci_wide()}, and \code{build_top_loci_export()} trio.
#' There is no separately exposed long or wide table.
#'
#' The output column order is exactly (22 columns):
#' \code{#chr}, \code{start}, \code{end}, \code{a1}, \code{a2},
#' \code{variant}, \code{gene}, \code{event},
#' \code{n}, \code{maf}, \code{beta}, \code{se},
#' \code{pip}, \code{posterior_effect_mean}, \code{posterior_effect_se},
#' \code{cs_95}, \code{cs_70}, \code{cs_50}, \code{cs_95_purity},
#' \code{method}, \code{grange_start}, \code{grange_end}.
#'
#' The \code{cs_95}, \code{cs_70}, \code{cs_50} columns are character strings of
#' the form \code{"<method>_<cs_index>"} where each method numbers its credible
#' sets independently starting at 1. Variants retained by the PIP cutoff but not
#' assigned to any credible set at the given coverage use \code{"<method>_0"}.
#' \code{cs_95_purity} is the 0.95-coverage credible-set purity for the row's
#' \code{(method, cs_95)}; rows whose \code{cs_95} is \code{"<method>_0"} carry
#' \code{cs_95_purity = 0}.
#'
#' Row uniqueness inside this function's output is one row per
#' \code{(variant, gene, cs_membership)} at the given \code{method}. Overlapping
#' credible-set membership for the same method produces one row per CS, so the
#' overlapping-CS contract is preserved.
#'
#' @param fit A fitted SuSiE-family object (must expose \code{alpha},
#'   \code{mu}, \code{mu2}, \code{pip}).
#' @param cs_tables A list of CS tables (one per coverage), as produced by
#'   \code{compute_cs_tables()}.
#' @param variant_names Character vector of variant IDs in
#'   \code{chr:pos:A2:A1} form, length equal to the number of variants in the
#'   fit. Used to construct \code{variant}, \code{#chr}, \code{start},
#'   \code{end}, \code{a1}, \code{a2}.
#' @param sumstats Optional marginal-association summary statistics
#'   (\code{betahat}, \code{sebetahat}) used to fill \code{beta} and \code{se}.
#' @param maf Optional numeric vector of minor-allele frequencies.
#' @param method Method name (e.g. \code{"susie"}, \code{"susie_inf"}). Used
#'   to construct the per-method \code{"<method>_<cs_index>"} strings and the
#'   \code{method} column. Required.
#' @param signal_cutoff PIP cutoff for retaining PIP-only (non-CS) variants.
#' @param data_x Optional regional genotype matrix; used only for sample-count
#'   shape checks.
#' @param data_y Optional regional phenotype matrix; \code{nrow(data_y)} fills
#'   \code{n}, \code{colnames(data_y)[1]} fills \code{gene}.
#' @param other_quantities Optional list with reserved subfields
#'   \code{region} (e.g. \code{"chr1:100-200"}) and \code{condition_id}
#'   (e.g. \code{"Ast_DeJager_eQTL"}); used to fill \code{grange_start},
#'   \code{grange_end}, and to compose \code{event} as
#'   \code{paste(condition_id, gene, sep = "_")}. Missing subfields are
#'   filled with \code{NA} rather than dropped.
#' @return A data frame in the fixed 22-column unified \code{top_loci} shape
#'   for this one fit and one method. Returns an empty data frame with the
#'   correct columns and dtypes if there is nothing to retain.
#' @export
build_top_loci <- function(fit, cs_tables, variant_names, sumstats = NULL,
                           maf = NULL, method, signal_cutoff = 0.1,
                           data_x = NULL, data_y = NULL,
                           other_quantities = NULL) {
  if (missing(method) || is.null(method) ||
      length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("build_top_loci: `method` is required (e.g. \"susie\", \"susie_inf\").")
  }
  if (length(cs_tables) == 0) return(.empty_top_loci())
  coverage_values <- attr(cs_tables, "coverage")
  if (is.null(coverage_values)) coverage_values <- rep(NA_real_, length(cs_tables))

  # Per-fit constants.
  data_y_mat <- if (!is.null(data_y)) as.matrix(data_y) else NULL
  fit_n    <- if (is.null(data_y_mat)) NA_integer_ else as.integer(nrow(data_y_mat))
  fit_gene <- if (!is.null(data_y_mat) && !is.null(colnames(data_y_mat))) {
    colnames(data_y_mat)[1]
  } else NA_character_
  fit_condition_id <- other_quantities$condition_id
  fit_event <- if (!is.null(fit_condition_id) &&
                   !is.na(fit_gene) && nzchar(fit_gene)) {
    paste(fit_condition_id, fit_gene, sep = "_")
  } else NA_character_
  grange_se        <- .parse_grange(other_quantities$region)
  fit_grange_start <- grange_se[["start"]]
  fit_grange_end   <- grange_se[["end"]]

  # Per-variant posterior effect and SE, computed once for all variants.
  alpha <- as.matrix(fit$alpha)
  mu    <- if (!is.null(fit$mu))  as.matrix(fit$mu)  else NULL
  mu2   <- if (!is.null(fit$mu2)) as.matrix(fit$mu2) else NULL
  posterior_effect <- if (!is.null(mu) && all(dim(alpha) == dim(mu))) {
    colSums(alpha * mu)
  } else rep(NA_real_, length(variant_names))
  posterior_effect_se <- if (!is.null(mu2) && all(dim(alpha) == dim(mu2))) {
    sqrt(pmax(colSums(alpha * mu2) - posterior_effect^2, 0))
  } else rep(NA_real_, length(variant_names))

  # Per-coverage credible-set purity vectors, indexed by 1-based CS index.
  purity_per_cov <- lapply(cs_tables, function(ct) {
    sets_purity <- ct$sets$purity
    if (!is.null(sets_purity) &&
        "min.abs.corr" %in% names(sets_purity)) {
      as.numeric(sets_purity$min.abs.corr)
    } else if (!is.null(ct$cs_corr)) {
      vapply(seq_along(ct$cs_corr), function(j) {
        m <- ct$cs_corr[[j]]
        if (is.null(m)) return(NA_real_)
        if (!is.matrix(m) || nrow(m) <= 1) return(1)
        min(abs(m[upper.tri(m)]))
      }, numeric(1))
    } else {
      rep(NA_real_, length(ct$sets$cs))
    }
  })

  # Internal long-shaped collection: one row per
  # (variant_idx, cs_idx_at_this_coverage, coverage). Not exposed. Used only
  # to project to the 22-column shape below.
  long_rows <- list()
  for (i in seq_along(cs_tables)) {
    ct <- cs_tables[[i]]
    top_variants_idx <- get_top_variants_idx(ct, signal_cutoff)
    cs_info <- get_cs_info(ct$sets$cs, top_variants_idx)
    if (is.null(cs_info) || nrow(cs_info) == 0) next
    long_rows[[length(long_rows) + 1L]] <- data.frame(
      variant_idx = as.integer(cs_info$variant_idx),
      cs_idx      = as.integer(cs_info$cs_idx),
      coverage    = as.numeric(coverage_values[[i]]),
      stringsAsFactors = FALSE
    )
  }
  if (length(long_rows) == 0) return(.empty_top_loci())
  long_df <- do.call(rbind, long_rows)
  if (nrow(long_df) == 0) return(.empty_top_loci())

  # Key grid: one row per (variant_idx, cs_idx) at the export grain. Preserves
  # overlapping CS membership for the same method.
  key_grid <- unique(long_df[, c("variant_idx", "cs_idx"), drop = FALSE])
  rownames(key_grid) <- NULL
  n_keys <- nrow(key_grid)

  # Helper: did this (variant_idx, cs_idx) appear at this coverage? If yes,
  # return cs_idx; otherwise 0.
  lookup_cs_at_cov <- function(v_idx, c_idx, cov) {
    sel <- long_df$variant_idx == v_idx &
           long_df$cs_idx       == c_idx &
           abs(long_df$coverage - cov) < 1e-12
    if (any(sel)) as.integer(c_idx) else 0L
  }
  cov95_table_idx <- which(abs(coverage_values - 0.95) < 1e-12)

  format_cs_string <- function(idx) {
    if (is.na(idx) || idx <= 0L) paste0(method, "_0") else paste0(method, "_", idx)
  }

  cs_95        <- character(n_keys)
  cs_70        <- character(n_keys)
  cs_50        <- character(n_keys)
  cs_95_purity <- numeric(n_keys)
  for (k in seq_len(n_keys)) {
    v_idx <- key_grid$variant_idx[k]
    c_idx <- key_grid$cs_idx[k]
    cs95_idx <- lookup_cs_at_cov(v_idx, c_idx, 0.95)
    cs70_idx <- lookup_cs_at_cov(v_idx, c_idx, 0.70)
    cs50_idx <- lookup_cs_at_cov(v_idx, c_idx, 0.50)
    cs_95[k] <- format_cs_string(cs95_idx)
    cs_70[k] <- format_cs_string(cs70_idx)
    cs_50[k] <- format_cs_string(cs50_idx)
    if (cs95_idx > 0L && length(cov95_table_idx) > 0L) {
      pvec <- purity_per_cov[[cov95_table_idx[1]]]
      cs_95_purity[k] <- if (cs95_idx <= length(pvec)) {
        val <- pvec[cs95_idx]
        if (is.na(val)) 0 else as.numeric(val)
      } else 0
    } else {
      cs_95_purity[k] <- 0
    }
  }

  # Per-variant lookups indexed by the row's variant_idx.
  v_idx_vec       <- key_grid$variant_idx
  variant_id_vec  <- variant_names[v_idx_vec]
  pip_vec         <- as.numeric(fit$pip[v_idx_vec])
  post_mean_vec   <- posterior_effect[v_idx_vec]
  post_se_vec     <- posterior_effect_se[v_idx_vec]
  beta_vec <- if (!is.null(sumstats$betahat))   sumstats$betahat[v_idx_vec]   else rep(NA_real_, n_keys)
  se_vec   <- if (!is.null(sumstats$sebetahat)) sumstats$sebetahat[v_idx_vec] else rep(NA_real_, n_keys)
  maf_vec  <- if (!is.null(maf))                maf[v_idx_vec]                else rep(NA_real_, n_keys)

  parsed <- tryCatch(
    suppressWarnings(parse_variant_id(variant_id_vec)),
    error = function(e) {
      stop("build_top_loci: parse_variant_id failed: ", conditionMessage(e))
    }
  )
  if (is.null(parsed) || nrow(parsed) != length(variant_id_vec)) {
    stop("build_top_loci: parse_variant_id did not return one row per variant.")
  }
  invalid <- is.na(parsed$chrom) | is.na(parsed$pos) |
    is.na(parsed$A1) | !nzchar(parsed$A1) |
    is.na(parsed$A2) | !nzchar(parsed$A2)
  if (any(invalid)) {
    first_bad <- variant_id_vec[which(invalid)[[1]]]
    stop("build_top_loci: parse_variant_id produced invalid coordinates ",
         "for variant_id: ", first_bad)
  }

  out <- data.frame(
    "#chr"                = parsed$chrom,
    start                 = as.integer(parsed$pos) - 1L,
    end                   = as.integer(parsed$pos),
    a1                    = parsed$A1,
    a2                    = parsed$A2,
    variant               = variant_id_vec,
    gene                  = rep(fit_gene, n_keys),
    event                 = rep(fit_event, n_keys),
    n                     = rep(fit_n, n_keys),
    maf                   = maf_vec,
    beta                  = beta_vec,
    se                    = se_vec,
    pip                   = pip_vec,
    posterior_effect_mean = post_mean_vec,
    posterior_effect_se   = post_se_vec,
    cs_95                 = cs_95,
    cs_70                 = cs_70,
    cs_50                 = cs_50,
    cs_95_purity          = cs_95_purity,
    method                = rep(method, n_keys),
    grange_start          = rep(fit_grange_start, n_keys),
    grange_end            = rep(fit_grange_end, n_keys),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
  rownames(out) <- NULL
  out
}

.empty_top_loci <- function() {
  data.frame(
    "#chr"                = character(),
    start                 = integer(),
    end                   = integer(),
    a1                    = character(),
    a2                    = character(),
    variant               = character(),
    gene                  = character(),
    event                 = character(),
    n                     = integer(),
    maf                   = numeric(),
    beta                  = numeric(),
    se                    = numeric(),
    pip                   = numeric(),
    posterior_effect_mean = numeric(),
    posterior_effect_se   = numeric(),
    cs_95                 = character(),
    cs_70                 = character(),
    cs_50                 = character(),
    cs_95_purity          = numeric(),
    method                = character(),
    grange_start          = integer(),
    grange_end            = integer(),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
}

.parse_grange <- function(region_str) {
  if (is.null(region_str) || length(region_str) == 0L ||
      is.na(region_str) || !nzchar(as.character(region_str))) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  pr <- tryCatch(parse_region(as.character(region_str)),
                 error = function(e) NULL)
  if (is.null(pr) || !is.data.frame(pr)) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  c(start = as.integer(pr$start), end = as.integer(pr$end))
}

# Project the new 22-column `top_loci` into the legacy shape expected by the
# FineMappingResult S4 slot, vcf_writer, getPIP, and getCS. We add backward-
# compatible aliases without renaming any column in the wrapper-facing
# `top_loci`:
#
#   * `variant_id` — copy of `variant`
#   * `cs`         — integer credible-set index derived from `cs_95` strings of
#                    the form `<method>_<idx>` (PIP-only `<method>_0` -> 0L)
#
# This isolates the schema change to susie_wrapper.R so AllClasses.R,
# AllMethods.R, and vcf_writer.R do not have to change.
.top_loci_for_s4_slot <- function(top_loci) {
  if (is.null(top_loci) || nrow(top_loci) == 0) {
    return(data.frame(variant_id = character(0),
                      method     = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- top_loci
  if ("variant" %in% names(out) && !"variant_id" %in% names(out)) {
    out$variant_id <- out$variant
  }
  if ("cs_95" %in% names(out) && !"cs" %in% names(out)) {
    out$cs <- vapply(out$cs_95, function(s) {
      if (is.na(s) || !nzchar(s)) return(0L)
      tail_str <- sub("^.*_", "", s)
      suppressWarnings(as.integer(tail_str))
    }, integer(1))
    out$cs[is.na(out$cs)] <- 0L
  }
  out
}

trim_finemapping_fit <- function(fit, effect_idx, method, cs_tables) {
  alpha <- .as_effect_matrix(fit$alpha)
  lbf_variable <- .as_lbf_matrix(fit)
  primary <- cs_tables[[1]]
  secondary <- if (length(cs_tables) > 1) {
    lapply(cs_tables[-1], function(x) x[names(x) != "pip"])
  } else {
    NULL
  }

  trimmed <- list(
    pip = as.numeric(fit$pip),
    sets = primary$sets,
    cs_corr = primary$cs_corr,
    sets_secondary = secondary,
    alpha = alpha[effect_idx, , drop = FALSE],
    lbf_variable = if (!is.null(lbf_variable)) lbf_variable[effect_idx, , drop = FALSE] else NULL,
    V = if (!is.null(fit$V)) fit$V[effect_idx] else NULL,
    niter = fit$niter,
    max_L = nrow(alpha),
    n_effects = nrow(alpha)
  )

  if (!is.null(fit$X_column_scale_factors)) trimmed$X_column_scale_factors <- fit$X_column_scale_factors
  if (!is.null(fit$mu)) {
    trimmed$mu <- if (length(dim(fit$mu)) == 3) fit$mu[effect_idx, , , drop = FALSE] else fit$mu[effect_idx, , drop = FALSE]
  }
  if (!is.null(fit$mu2)) {
    # mu2 is L x p for univariate susie and L x p x R for multivariate (mvsusie).
    # Match the shape handling used for mu just above.
    trimmed$mu2 <- if (length(dim(fit$mu2)) == 3) fit$mu2[effect_idx, , , drop = FALSE] else fit$mu2[effect_idx, , drop = FALSE]
  }
  if (!is.null(fit$theta)) trimmed$theta <- fit$theta
  if (!is.null(fit$omega_weights)) trimmed$omega_weights <- fit$omega_weights

  if (method == "mvsusie") {
    if (!is.null(fit$mu2_diag)) trimmed$mu2_diag <- fit$mu2_diag[effect_idx, , , drop = FALSE]
    if (requireNamespace("mvsusieR", quietly = TRUE)) {
      trimmed$coef <- mvsusieR::coef.mvsusie(fit)[-1, , drop = FALSE]
    }
    if (!is.null(fit$conditional_lfsr)) trimmed$clfsr <- fit$conditional_lfsr[effect_idx, , , drop = FALSE]
  }

  class(trimmed) <- unique(c(method, "susie"))
  trimmed
}

#' Format Fine-mapping Post-processing for Protocol Output
#'
#' Converts method-aware fine-mapping post-processing output into the root-level
#' fields consumed by protocol RDS files. The single top-loci output is the
#' \code{top_loci} field (the 22-column unified table); there is no
#' \code{top_loci_long}, no wide-format \code{top_loci}, and no
#' \code{top_loci_export}.
#'
#' @param post Output from \code{\link{postprocess_finemapping_fits}}.
#' @param primary_method Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{variant_names},
#'   \code{susie_result_trimmed}, and the single unified \code{top_loci}
#'   22-column table.
#' @export
format_finemapping_output <- function(post, primary_method) {
  method_post <- post$finemapping_results[[primary_method]]
  if (is.null(method_post)) {
    stop("primary_method was not found in finemapping_results: ", primary_method)
  }
  keep_names <- setdiff(names(method_post), c("result_trimmed", "top_loci"))
  c(
    method_post[keep_names],
    list(
      susie_result_trimmed = method_post$result_trimmed,
      top_loci = post$top_loci
    )
  )
}

#' Adjust SuSiE Weights
#'
#' Adjusts SuSiE TWAS weights by subsetting to intersected variants and
#' optionally running allele QC against LD reference variants.
#'
#' @param twas_weights_results A list containing TWAS weight data (nested structure).
#' @param keep_variants Vector of variant names to keep.
#' @param run_allele_qc Whether to run allele_qc to align alleles. Default TRUE.
#' @param variable_name_obj Path to variant names in the nested list.
#' @param susie_obj Path to susie result in the nested list.
#' @param twas_weights_table Path to weights table in the nested list.
#' @param LD_variants Vector of LD reference variant IDs for allele QC.
#' @param match_min_prop Minimum proportion of matched variants. Default 0.2.
#' @return A list with adjusted_susie_weights and remained_variants_ids.
#' @export
adjust_susie_weights <- function(twas_weights_results, keep_variants, run_allele_qc = TRUE,
                                 variable_name_obj = c("susie_results", context, "variant_names"),
                                 susie_obj = c("susie_results", context, "susie_result_trimmed"),
                                 twas_weights_table = c("weights", context), LD_variants, match_min_prop = 0.2) {
  # Intersect the rownames of weights with keep_variants
  twas_weights_variants <- get_nested_element(twas_weights_results, variable_name_obj)
  # Normalize to canonical format (with chr prefix)
  twas_weights_variants <- normalize_variant_id(twas_weights_variants)
  # allele flip twas weights matrix variants name
  if (run_allele_qc) {
    weights_matrix <- get_nested_element(twas_weights_results, twas_weights_table)
    if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weights_matrix))) {
      weights_matrix <- cbind(parse_variant_id(twas_weights_variants), weights_matrix)
    }
    weights_matrix_qced <- match_ref_panel(weights_matrix, LD_variants, colnames(weights_matrix)[!colnames(weights_matrix) %in% c(
      "chrom",
      "pos", "A2", "A1"
    )], match_min_prop = match_min_prop)
    # match_ref_panel outputs canonical variant_ids (with chr prefix)
    original_idx <- match(weights_matrix_qced$qc_summary$variants_id_original, twas_weights_variants)
    intersected_indices <- original_idx[weights_matrix_qced$qc_summary$keep == TRUE]
  } else {
    # Normalize keep_variants to canonical format for matching
    keep_variants_normalized <- normalize_variant_id(keep_variants)
    intersected_variants <- intersect(twas_weights_variants, keep_variants_normalized)
    intersected_indices <- match(intersected_variants, twas_weights_variants)
  }
  if (length(intersected_indices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbf_matrix <- get_nested_element(twas_weights_results, c(susie_obj, "lbf_variable"))
  mu <- get_nested_element(twas_weights_results, c(susie_obj, "mu"))
  x_column_scal_factors <- get_nested_element(twas_weights_results, c(susie_obj, "X_column_scale_factors"))

  lbf_matrix_subset <- lbf_matrix[, intersected_indices, drop = FALSE]
  mu_subset <- mu[, intersected_indices, drop = FALSE]
  x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
  adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset
  # allele_qc now outputs canonical variant_ids (with chr prefix) -- no need to add chr
  remained_variants_ids <- if (run_allele_qc) {
    weights_matrix_qced$target_data_qced$variant_id
  } else {
    intersected_variants
  }
  return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = remained_variants_ids))
}

#' Run the SuSiE RSS pipeline
#'
#' Runs SuSiE RSS analysis with the specified method. Supports both z+R
#' (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstats Data frame with 'z' or ('beta' and 'se') columns.
#' @param LD_mat LD correlation matrix. Mutually exclusive with X_mat.
#' @param X_mat Genotype matrix (samples x variants). Mutually exclusive with LD_mat.
#' @param n Sample size.
#' @param L Maximum number of causal configurations (default: 30).
#' @param L_greedy Initial greedy number of causal configurations (default: 5).
#' @param analysis_method One of "susie_rss", "single_effect", "bayesian_conditional_regression".
#' @param coverage Coverage level (default: 0.95).
#' @param secondary_coverage Secondary coverage levels (default: c(0.7, 0.5)).
#' @param signal_cutoff PIP cutoff for selecting top loci (default: 0.1).
#' @param min_abs_corr Minimum absolute correlation for CS purity (default: 0.8).
#' @param R_finite Controls variance inflation to account for estimating
#'   the R matrix from a finite reference panel. NULL (default): no
#'   variance inflation. Passed directly to susie_rss.
#' @param R_mismatch LD mismatch correction method passed directly to susie_rss.
#'   Default NULL disables mismatch correction.
#' @param ... Additional parameters passed to susie_rss (e.g., var_y).
#' @return A list with post-processed SuSiE RSS results. Method-specific
#'   \code{top_loci} columns use the selected \code{analysis_method}, for
#'   example \code{pip_susie_rss} or \code{pip_single_effect}.
#' @importFrom susieR susie_rss
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susie_rss_pipeline <- function(sumstats, LD_mat = NULL, X_mat = NULL, n = NULL,
                               L = 30, L_greedy = 5,
                               analysis_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
                               coverage = 0.95,
                               secondary_coverage = c(0.7, 0.5),
                               signal_cutoff = 0.1,
                               min_abs_corr = 0.8,
                               R_finite = NULL, R_mismatch = NULL, ...) {
  analysis_method <- match.arg(analysis_method)
  if (is.null(LD_mat) && is.null(X_mat)) stop("Either LD_mat or X_mat must be provided.")
  if (!is.null(LD_mat) && !is.null(X_mat)) stop("Only one of LD_mat or X_mat should be provided, not both.")
  if (!is.null(L_greedy)) L_greedy <- min(L_greedy, L)

  if (!is.null(sumstats$z)) {
    z <- sumstats$z
  } else if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
    z <- sumstats$beta / sumstats$se
  } else {
    stop("sumstats must have 'z' or ('beta' and 'se') columns.")
  }
  if (is.null(names(z)) && !is.null(sumstats$variant_id) && length(sumstats$variant_id) == length(z)) {
    names(z) <- sumstats$variant_id
  }
  if (is.null(names(z)) && !is.null(rownames(sumstats)) && length(rownames(sumstats)) == length(z)) {
    names(z) <- rownames(sumstats)
  }

  common <- list(z = z, n = n, coverage = coverage,
                 R_finite = R_finite, R_mismatch = R_mismatch, ...)
  if (!is.null(X_mat)) common$X <- X_mat else common$R <- LD_mat

  if (analysis_method == "single_effect") {
    res <- do.call(susie_rss, c(common, list(L = 1, L_greedy = NULL, max_iter = 1)))
  } else if (analysis_method == "bayesian_conditional_regression") {
    res <- do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy, max_iter = 1)))
  } else {
    res <- do.call(susie_rss, c(common, list(L = L, L_greedy = L_greedy)))
  }

  # For post-processing, need a square matrix (R or computed from X).
  # For mixture panels (list of X), use the first panel to compute R.
  if (!is.null(LD_mat)) {
    data_x <- LD_mat
  } else if (is.list(X_mat) && !is.matrix(X_mat)) {
    data_x <- compute_LD(X_mat[[1]][, seq_along(z), drop = FALSE], method = "sample")
  } else {
    data_x <- compute_LD(X_mat[, seq_along(z), drop = FALSE], method = "sample")
  }

  rss_method <- analysis_method
  rss_fit <- .set_finemapping_fit_class(res, rss_method)
  post <- postprocess_finemapping_fits(
    fits = setNames(list(rss_fit), rss_method),
    data_x = data_x,
    data_y = list(z = z),
    coverage = coverage,
    secondary_coverage = secondary_coverage,
    signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr
  )
  format_finemapping_output(post, primary_method = rss_method)
}

#' @noRd
get_cs_index <- function(snps_idx, susie_cs) {
  # Return ALL CS indices that contain this variant (not just one)
  idx <- which(vapply(susie_cs, function(x) snps_idx %in% x, logical(1)))
  if (length(idx) == 0) return(NA_integer_)
  return(idx)
}
#' @noRd
get_top_variants_idx <- function(susie_output, signal_cutoff) {
  c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>%
    unique() %>%
    sort()
}
# Returns a data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair.
# Variants in multiple CSs get multiple rows.
#' @importFrom stringr str_replace
#' @noRd
get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
  cs_names <- names(susie_output_sets_cs)
  rows <- lapply(top_variants_idx, function(vi) {
    idx <- get_cs_index(vi, susie_output_sets_cs)
    if (length(idx) == 1 && is.na(idx)) {
      data.frame(variant_idx = vi, cs_idx = 0L, stringsAsFactors = FALSE)
    } else {
      cs_nums <- as.integer(str_replace(cs_names[idx], "L", ""))
      data.frame(variant_idx = rep(vi, length(cs_nums)), cs_idx = cs_nums, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
