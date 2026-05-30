context("susie_finemapping")

# =============================================================================
# lbf_to_alpha_vector (internal)
# =============================================================================

test_that("lbf_to_alpha_vector converts correctly", {
  lbf <- c(a = -0.5, b = 1.2, c = 0.3)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 3)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with prior weights", {
  lbf <- c(a = 1, b = 1, c = 1)  # Equal LBFs
  pw <- c(0.5, 0.25, 0.25)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_true(alpha[1] > alpha[2])
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha_vector returns zeros for all-zero lbf", {
  lbf <- c(a = 0, b = 0, c = 0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(alpha == 0))
})

test_that("lbf_to_alpha_vector handles single element", {
  lbf <- c(a = 2.0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 1)
  expect_equal(alpha[["a"]], 1.0)
})

test_that("lbf_to_alpha_vector handles very large LBFs without overflow", {
  lbf <- c(a = 500, b = 500.1, c = 499)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector handles very negative LBFs", {
  lbf <- c(a = -1000, b = -999, c = -1001)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with unequal prior weights", {
  lbf <- c(a = 0.5, b = 0.5, c = 0.5)
  pw <- c(0.8, 0.1, 0.1)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha[1] > 0.7)
})

# =============================================================================
# lbf_to_alpha (matrix version)
# =============================================================================

test_that("lbf_to_alpha converts log BFs to posteriors", {
  lbf <- matrix(c(0, 3, 2, 1, 4, 0), nrow = 2, ncol = 3)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(dim(alpha), c(2, 3))
  expect_equal(rowSums(alpha), c(1, 1), tolerance = 1e-10)
  expect_true(alpha[1, 3] > alpha[1, 1])
  expect_true(alpha[2, 1] > alpha[2, 3])
})

test_that("lbf_to_alpha handles uniform lbf", {
  lbf <- matrix(1, nrow = 1, ncol = 5)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(as.numeric(alpha), rep(0.2, 5), tolerance = 1e-10)
})

test_that("lbf_to_alpha handles single-row matrix", {
  lbf <- matrix(c(1.0, 2.0, 0.5), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha handles large matrix", {
  set.seed(42)
  lbf <- matrix(rnorm(100), nrow = 10, ncol = 10)
  colnames(lbf) <- paste0("v", 1:10)
  result <- lbf_to_alpha(lbf)
  expect_equal(dim(result), c(10, 10))
  expect_equal(rowSums(result), rep(1, 10), tolerance = 1e-10)
})

test_that("lbf_to_alpha with mixed zero and nonzero rows", {
  lbf <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 2, byrow = TRUE)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_true(all(result[1, ] == 0))
  expect_equal(sum(result[2, ]), 1, tolerance = 1e-10)
})

# =============================================================================
# get_cs_index (internal)
# =============================================================================

test_that("get_cs_index finds variant in credible set", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(idx), 1)
})

test_that("get_cs_index returns NA for variant not in any CS", {
  susie_cs <- list(L1 = c(1, 2), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(3, susie_cs)
  expect_true(is.na(idx))
})

test_that("get_cs_index returns all CS indices when variant in multiple", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(2, 4, 5))
  idx <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(idx), c(1, 2))
})

test_that("get_cs_index returns all matching CS regardless of size", {
  susie_cs <- list(L1 = c(1, 2, 3, 4, 5), L2 = c(2, 3))
  result <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(result), c(1, 2))
})

test_that("get_cs_index handles empty CS list", {
  susie_cs <- list()
  result <- pecotmr:::get_cs_index(1, susie_cs)
  expect_true(is.na(result))
})

test_that("get_cs_index returns correct CS assignment with real susie fit", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  y <- X %*% beta + rnorm(n, sd = 0.5)
  fit <- susieR::susie(X, y, L = 5)
  # With beta[1]=2 and sd=0.5, susie should find a CS containing variant 1
  expect_false(is.null(fit$sets$cs))
  idx <- pecotmr:::get_cs_index(1, fit$sets$cs)
  expect_true(is.numeric(unname(idx)))
  expect_true(all(idx >= 1))
})

# =============================================================================
# get_top_variants_idx (internal)
# =============================================================================

test_that("get_top_variants_idx returns combined PIP and CS variants", {
  susie_output <- list(
    pip = c(0.01, 0.15, 0.02, 0.5, 0.01),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_true(1 %in% result)
  expect_true(2 %in% result)
  expect_true(4 %in% result)
  expect_true(all(result == sort(result)))
})

test_that("get_top_variants_idx with no CS", {
  susie_output <- list(
    pip = c(0.01, 0.5, 0.02, 0.8, 0.01),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_equal(result, c(2, 4))
})

test_that("get_top_variants_idx with all low PIPs", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_equal(result, c(1, 2))
})

test_that("get_top_variants_idx with high cutoff and no CS", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_length(result, 0)
})

# =============================================================================
# get_cs_info (internal)
# =============================================================================

test_that("get_cs_info maps variants to CS numbers", {
  susie_cs <- list(L1 = c(1, 2), L3 = c(4, 5, 6))
  top_idx <- c(1, 3, 5)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  # Now returns data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair
  expect_true(is.data.frame(result))
  expect_equal(result$variant_idx, c(1, 3, 5))
  expect_equal(result$cs_idx, c(1L, 0L, 3L))
})

test_that("get_cs_info handles all variants outside CS", {
  susie_cs <- list(L1 = c(1, 2))
  top_idx <- c(5, 6, 7)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  expect_true(all(result$cs_idx == 0))
})

test_that("get_cs_info reports variant in multiple CSs as multiple rows", {
  susie_cs <- list(L1 = c(1, 2, 3), L3 = c(2, 3, 4))
  top_idx <- c(1, 2, 4)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  # variant 2 is in both L1 and L3, so it gets two rows
  expect_equal(nrow(result), 4)
  expect_equal(sum(result$variant_idx == 2), 2)
  expect_equal(sort(result$cs_idx[result$variant_idx == 2]), c(1L, 3L))
})

# =============================================================================
# susie_rss_pipeline
# =============================================================================

test_that("susie_rss_pipeline errors on missing z and beta/se", {
  sumstats <- data.frame(x = 1)
  LD_mat <- matrix(1)
  expect_error(susie_rss_pipeline(sumstats, LD_mat), "must have 'z'")
})

test_that("susie_rss_pipeline errors on invalid method", {
  sumstats <- list(z = rnorm(5))
  LD_mat <- diag(5)
  expect_error(susie_rss_pipeline(sumstats, LD_mat, analysis_method = "invalid"))
})

test_that("susie_rss_pipeline runs with single_effect method", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R, analysis_method = "single_effect")
  expect_true(is.list(result))
  expect_true("variant_names" %in% names(result))
  expect_true("susie_result_trimmed" %in% names(result))
  if (!is.null(result$top_loci)) {
    expect_true("pip_single_effect" %in% names(result$top_loci))
    expect_true("CS_95_single_effect" %in% names(result$top_loci))
  }
  # PIPs should be numeric, in [0,1], and sum to at most 1 (L=1)
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 1 + 1e-6)
  # Credible sets, if any, should contain valid indices
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline runs with bayesian_conditional_regression", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "bayesian_conditional_regression",
    L = 5, L_greedy = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  if (!is.null(result$top_loci)) {
    expect_true("pip_bayesian_conditional_regression" %in% names(result$top_loci))
    expect_true("CS_95_bayesian_conditional_regression" %in% names(result$top_loci))
  }
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  # With L=5, sum of PIPs can be up to L
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline uses beta/se when z not provided", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 15
  beta <- rnorm(n, sd = 0.1)
  se <- rep(0.1, n)
  names(beta) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(beta)
  sumstats <- list(beta = beta, se = se)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "susie_rss",
    L = 5, L_greedy = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

# =============================================================================
# susie_weights
# =============================================================================

test_that("susie_weights returns zeros when fit lacks alpha/mu", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_weights(susie_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

test_that("susie_weights checks dimension mismatch", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  fake_fit <- list(pip = rep(0.01, 10))
  expect_error(susie_weights(X = X, susie_fit = fake_fit), "Dimension mismatch")
})

# =============================================================================
# susie_ash_weights
# =============================================================================

test_that("susie_ash_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_ash_weights(susie_ash_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# susie_inf_weights
# =============================================================================

test_that("susie_inf_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_inf_weights(susie_inf_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# glmnet_weights
# =============================================================================

test_that("glmnet_weights produces non-zero weights for correlated data", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, -2, rep(0, p - 2))
  y <- X %*% beta_true + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 0.5)
  expect_length(w, p)
  expect_true(any(w != 0))
})

test_that("glmnet_weights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 3] <- 1  # zero variance column
  y <- X[, 1] * 2 + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 1)
  expect_length(w, p)
  expect_equal(w[3], 0)
})

# =============================================================================
# init_prior_sd
# =============================================================================

test_that("init_prior_sd returns n standard deviations", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n_samples <- 50
  p <- 10
  X <- matrix(rnorm(n_samples * p), n_samples, p)
  y <- X[, 1] * 2 + rnorm(n_samples)

  sds <- pecotmr:::init_prior_sd(X, y, n = 15)
  expect_length(sds, 15)
  expect_equal(sds[1], 0)
  expect_true(all(diff(sds) >= 0))
})

# =============================================================================
# adjust_susie_weights
# =============================================================================

# Helper: build a minimal twas_weights_results object with the nested structure
# adjust_susie_weights expects (susie_results / weights paths).
make_adjust_obj <- function(variant_ids, L = 3, ctx = "ctx") {
  set.seed(123)
  p <- length(variant_ids)
  weights_df <- data.frame(
    susie = rnorm(p), enet = rnorm(p),
    row.names = variant_ids, stringsAsFactors = FALSE
  )
  list(
    susie_results = setNames(list(list(
      variant_names = variant_ids,
      susie_result_trimmed = list(
        lbf_variable = matrix(rnorm(L * p), nrow = L, ncol = p),
        mu = matrix(rnorm(L * p), nrow = L, ncol = p),
        X_column_scale_factors = rep(1, p)
      )
    )), ctx),
    weights = setNames(list(weights_df), ctx)
  )
}

# Use non-strand-ambiguous alleles (A2="A", A1="G") so allele_qc keeps them.
adjust_vids <- function(positions = 1:6) {
  paste0("chr1:", positions, ":A:G")
}

# =============================================================================
# postprocess_finemapping_fits: analysis_script and V=NULL branches (Tier 1)
# =============================================================================

# Helper: build a minimal synthetic SuSiE-family output for post-processing
make_fake_susie_output <- function(p = 5, L = 3, has_V = TRUE) {
  vnames <- paste0("chr1:", 1:p, ":A:G")
  out <- list(
    pip = setNames(rep(0.01, p), vnames),
    alpha = matrix(1 / p, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    sets = list(
      cs = NULL,
      requested_coverage = 0.95
    ),
    niter = 10
  )
  if (has_V) {
    out$V <- rep(1, L)
  }
  out
}

test_that("postprocess_finemapping_fits stores analysis_script when load_script returns non-empty", {
  skip_if_not_installed("susieR")
  p <- 5
  fake_output <- make_fake_susie_output(p)
  R <- diag(p)
  colnames(R) <- rownames(R) <- names(fake_output$pip)
  local_mocked_bindings(
    load_script = function() "fake_script_content"
  )
  post <- postprocess_finemapping_fits(
    fits = list(susie_rss = pecotmr:::.set_finemapping_fit_class(fake_output, "susie_rss")),
    data_x = R,
    data_y = list(z = rnorm(p)),
    coverage = 0.95
  )
  result <- format_finemapping_output(post, primary_method = "susie_rss")
  expect_equal(result$analysis_script, "fake_script_content")
})

test_that("postprocess_finemapping_fits keeps all effects when V is NULL", {
  skip_if_not_installed("susieR")
  p <- 5
  L <- 3
  fake_output <- make_fake_susie_output(p, L = L, has_V = FALSE)
  R <- diag(p)
  colnames(R) <- rownames(R) <- names(fake_output$pip)
  post <- postprocess_finemapping_fits(
    fits = list(susie_rss = pecotmr:::.set_finemapping_fit_class(fake_output, "susie_rss")),
    data_x = R,
    data_y = list(z = rnorm(p)),
    coverage = 0.95
  )
  result <- format_finemapping_output(post, primary_method = "susie_rss")
  # With V=NULL, eff_idx = 1:L, so trimmed alpha should keep all L rows
  expect_equal(nrow(result$susie_result_trimmed$alpha), L)
  # V should be NULL in trimmed output
  expect_null(result$susie_result_trimmed$V)
})

# =============================================================================
# postprocess_finemapping_fits: mvsusie output (outcome_names, coef, clfsr)
# =============================================================================

test_that("postprocess_finemapping_fits stores outcome_names, coef, and clfsr for mvsusie", {
  skip_if_not_installed("susieR")
  skip_if_not_installed("mvsusieR")
  p <- 5
  L <- 3
  R <- 2
  vnames <- paste0("chr1:", 1:p, ":A:G")
  cnames <- paste0("cond_", 1:R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  fake_output <- list(
    pip = setNames(rep(0.01, p), vnames),
    alpha = matrix(1 / p, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    sets = list(cs = NULL, requested_coverage = 0.95),
    niter = 10,
    V = rep(1, L),
    outcome_names = cnames,
    conditional_lfsr = array(0.5, dim = c(L, p, R))
  )

  n <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vnames

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  post <- postprocess_finemapping_fits(
    fits = list(mvsusie = pecotmr:::.set_finemapping_fit_class(fake_output, "mvsusie")),
    data_x = X,
    data_y = NULL,
    X_scalar = 1, y_scalar = 1,
    coverage = 0.95
  )
  result <- format_finemapping_output(post, primary_method = "mvsusie")

  # outcome_names should be stored as context_names
  expect_equal(result$context_names, cnames)
  # coef should come from mvsusieR::coef.mvsusie
  expect_equal(result$susie_result_trimmed$coef, fake_coef[-1, , drop = FALSE])
  # conditional_lfsr should be trimmed to eff_idx
  expect_equal(dim(result$susie_result_trimmed$clfsr), c(L, p, R))
})

# =============================================================================
# susie_rss_pipeline X-mode branches (Tier 2)
# =============================================================================

test_that("susie_rss_pipeline X-mode passes X to susie_rss and computes LD from X for post-processor", {
  skip_if_not_installed("susieR")
  p <- 5
  n <- 20
  z <- rnorm(p)
  vnames <- paste0("chr1:", 1:p, ":A:G")
  names(z) <- vnames
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vnames

  captured_susie_args <- NULL
  captured_pp_data_x <- NULL
  local_mocked_bindings(
    susie_rss = function(...) {
      captured_susie_args <<- list(...)
      list(
        pip = setNames(rep(0.01, p), vnames),
        alpha = matrix(1 / p, nrow = 5, ncol = p),
        lbf_variable = matrix(0, nrow = 5, ncol = p),
        V = rep(1, 5),
        sets = list(cs = NULL, requested_coverage = 0.95),
        niter = 10
      )
    },
    postprocess_finemapping_fits = function(fits, data_x, ...) {
      captured_pp_data_x <<- data_x
      list()
    },
    format_finemapping_output = function(post, primary_method) list(variant_names = vnames)
  )
  result <- susie_rss_pipeline(list(z = z), X_mat = X, R_mismatch = "eb")
  expect_true("X" %in% names(captured_susie_args))
  expect_null(captured_susie_args$R)
  expect_equal(captured_susie_args$R_mismatch, "eb")
  # Post-processor should have received a p x p matrix (LD computed from X)
  expect_equal(dim(captured_pp_data_x), c(p, p))
})

# =============================================================================
# susie_rss_pipeline: mixture-panel (list of X) branch (Tier 2)
# =============================================================================

test_that("susie_rss_pipeline computes LD from first panel when X_mat is a list", {
  skip_if_not_installed("susieR")
  p <- 5
  n1 <- 20
  n2 <- 15
  z <- rnorm(p)
  vnames <- paste0("chr1:", 1:p, ":A:G")
  names(z) <- vnames
  X1 <- matrix(rnorm(n1 * p), n1, p)
  X2 <- matrix(rnorm(n2 * p), n2, p)
  colnames(X1) <- colnames(X2) <- vnames
  X_list <- list(panel1 = X1, panel2 = X2)

  captured_pp_data_x <- NULL
  local_mocked_bindings(
    susie_rss = function(...) {
      list(
        pip = setNames(rep(0.01, p), vnames),
        alpha = matrix(1 / p, nrow = 5, ncol = p),
        lbf_variable = matrix(0, nrow = 5, ncol = p),
        V = rep(1, 5),
        sets = list(cs = NULL, requested_coverage = 0.95),
        niter = 10
      )
    },
    postprocess_finemapping_fits = function(fits, data_x, ...) {
      captured_pp_data_x <<- data_x
      list()
    },
    format_finemapping_output = function(post, primary_method) list(variant_names = vnames)
  )
  result <- susie_rss_pipeline(list(z = z), X_mat = X_list)
  # data_x should be a p x p correlation matrix computed from X1 (first panel)
  expect_equal(dim(captured_pp_data_x), c(p, p))
  # It should be a symmetric matrix (correlation/LD)
  expect_equal(captured_pp_data_x, t(captured_pp_data_x), tolerance = 1e-10)
})

# =============================================================================
# adjust_susie_weights
# =============================================================================

test_that("adjust_susie_weights errors when no variants intersect", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  expect_error(
    adjust_susie_weights(
      obj,
      keep_variants = paste0("chr2:", 1:5, ":A:G"),
      run_allele_qc = FALSE,
      variable_name_obj = c("susie_results", "ctx", "variant_names"),
      susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
      twas_weights_table = c("weights", "ctx"),
      LD_variants = NULL
    ),
    "No intersected variants"
  )
})

test_that("adjust_susie_weights run_allele_qc=FALSE returns intersect coefs", {
  vids <- adjust_vids(1:6)
  obj <- make_adjust_obj(vids)
  keep <- vids[2:5]
  out <- adjust_susie_weights(
    obj,
    keep_variants = keep, run_allele_qc = FALSE,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    LD_variants = NULL
  )
  expect_length(out$adjusted_susie_weights, 4)
  expect_equal(out$remained_variants_ids, normalize_variant_id(keep))
  expect_true(all(is.finite(out$adjusted_susie_weights)))
})

test_that("adjust_susie_weights run_allele_qc=FALSE normalizes variant ids before matching", {
  # Object has non-canonical (no chr prefix) variant ids
  vids_raw <- c("1:1:A:G", "1:2:A:G", "1:3:A:G", "1:4:A:G")
  obj <- make_adjust_obj(vids_raw)
  # keep_variants supplied with chr prefix
  keep <- c("chr1:2:A:G", "chr1:3:A:G")
  out <- adjust_susie_weights(
    obj,
    keep_variants = keep, run_allele_qc = FALSE,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    LD_variants = NULL
  )
  expect_length(out$adjusted_susie_weights, 2)
  expect_equal(out$remained_variants_ids, c("chr1:2:A:G", "chr1:3:A:G"))
})

test_that("adjust_susie_weights run_allele_qc=TRUE returns adjusted xQTL coefs", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  out <- adjust_susie_weights(
    obj,
    keep_variants = vids, run_allele_qc = TRUE, LD_variants = vids,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    match_min_prop = 0.1
  )
  expect_true(length(out$adjusted_susie_weights) > 0)
  expect_true(all(is.finite(out$adjusted_susie_weights)))
  expect_true(all(grepl("^chr1:", out$remained_variants_ids)))
})

test_that("adjust_susie_weights run_allele_qc=TRUE auto-prepends chrom/pos/A2/A1", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  # Confirm the helper produced a weights matrix WITHOUT chrom/pos/A2/A1 cols
  expect_false(any(c("chrom", "pos", "A2", "A1") %in% colnames(obj$weights$ctx)))
  out <- adjust_susie_weights(
    obj,
    keep_variants = vids, run_allele_qc = TRUE, LD_variants = vids,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    match_min_prop = 0.1
  )
  expect_true(length(out$adjusted_susie_weights) > 0)
})

test_that("format_finemapping_output does not duplicate top loci variants", {
  top_loci <- data.frame(
    variant_id = paste0("v", 1:4),
    CS_95_susie = c(0L, 1L, NA_integer_, 0L),
    pip_susie = c(0.2, 0.005, 0.001, 0),
    stringsAsFactors = FALSE
  )
  post <- list(
    finemapping_results = list(susie = list(
      variant_names = paste0("v", 1:4),
      result_trimmed = list(pip = 1:4),
      top_loci_long = NULL
    )),
    top_loci_long = NULL,
    top_loci = top_loci
  )
  out <- format_finemapping_output(post, "susie")
  expect_false("top_loci_variants" %in% names(out))
  expect_equal(unique(out$top_loci$variant_id), paste0("v", 1:4))
})

.make_univariate_data <- function(seed = 42, n = 300, p = 50,
                                  effect_idx = integer(0), effect_size = NULL) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- sprintf("chr1:%d:G:A", seq_len(p))
  beta <- rep(0, p)
  if (length(effect_idx) > 0) {
    if (is.null(effect_size)) effect_size <- rep(1.5, length(effect_idx))
    beta[effect_idx] <- effect_size
  }
  y <- as.numeric(X %*% beta) + rnorm(n, sd = 0.5)
  list(X = X, y = y)
}

test_that("postprocess top_loci: single susie with signal yields per-method columns and no unsuffixed pip", {
  d <- .make_univariate_data(seed = 11, effect_idx = c(10, 30))
  fit <- susieR::susie(d$X, d$y, L = 5)
  post <- postprocess_finemapping_fits(list(susie = fit), data_x = d$X, data_y = d$y, coverage = 0.95)
  expect_false("pip" %in% colnames(post$top_loci))
  expect_true("pip_susie" %in% colnames(post$top_loci))
  expect_true("CS_95_susie" %in% colnames(post$top_loci))
  # No other method's columns
  expect_length(grep("_susie_inf$", colnames(post$top_loci)), 0)
  expect_length(grep("_mvsusie$",  colnames(post$top_loci)), 0)
})

test_that("postprocess top_loci: susie + susie_inf both with signal yield symmetric pip/CS columns", {
  d <- .make_univariate_data(seed = 12, effect_idx = c(15, 45))
  fits <- fit_susie_inf_then_susie(d$X, d$y)
  expect_gt(max(fits$susie$pip), 0.5)
  expect_gt(max(fits$susie_inf$pip), 0.5)
  post <- postprocess_finemapping_fits(fits, data_x = d$X, data_y = d$y, coverage = 0.95)
  cols <- colnames(post$top_loci)
  expect_false("pip" %in% cols)
  for (m in c("susie", "susie_inf")) {
    expect_true(paste0("pip_", m)        %in% cols, info = m)
    expect_true(paste0("CS_95_", m)      %in% cols, info = m)
    expect_true(paste0("CS_70_", m)      %in% cols, info = m)
    expect_true(paste0("CS_50_", m)      %in% cols, info = m)
  }
  # top_loci_long contains both methods
  expect_setequal(unique(post$top_loci_long$method), c("susie", "susie_inf"))
})

test_that("postprocess top_loci: when one method has no rows in long, its pip column still appears", {
  # Construct a case where susie has signal but the susie_inf trim gives no CS rows.
  # We use a high signal_cutoff so susie_inf may not produce top variants.
  d <- .make_univariate_data(seed = 13, effect_idx = c(20))
  fits <- fit_susie_inf_then_susie(d$X, d$y)
  post <- postprocess_finemapping_fits(fits, data_x = d$X, data_y = d$y,
                                       coverage = 0.95, signal_cutoff = 0.99)
  cols <- colnames(post$top_loci)
  # Both pip columns are present because they come from result_trimmed$pip,
  # independent of whether the method contributed rows to top_loci_long.
  expect_true("pip_susie" %in% cols)
  expect_true("pip_susie_inf" %in% cols)
  expect_false("pip" %in% cols)
})

test_that("postprocess top_loci: susie + susie_inf + mvsusie all run through trim without dimension errors", {
  skip_if_not_installed("mvsusieR")
  d <- .make_univariate_data(seed = 14, effect_idx = c(15, 45))
  Y_mv <- cbind(d$y, as.numeric(d$X %*% rep(0, ncol(d$X))) + rnorm(length(d$y), sd = 0.5))
  colnames(Y_mv) <- c("ctx1", "ctx2")
  fits_uni <- fit_susie_inf_then_susie(d$X, d$y)
  mv_fit <- mvsusieR::mvsusie(d$X, Y_mv, L = 5,
                              prior_variance = mvsusieR::create_mixture_prior(R = ncol(Y_mv)),
                              prior_weights = rep(1 / ncol(d$X), ncol(d$X)))
  fits_all <- list(susie = fits_uni$susie,
                   susie_inf = fits_uni$susie_inf,
                   mvsusie = mv_fit)
  # Bug B regression: this used to error with "incorrect number of dimensions"
  # because trim_finemapping_fit indexed mu2 as 2-D regardless of mvsusie's 3-D.
  expect_no_error(
    post <- postprocess_finemapping_fits(fits_all, data_x = d$X, data_y = Y_mv, coverage = 0.95)
  )
  for (m in c("susie", "susie_inf", "mvsusie")) {
    expect_true(paste0("pip_", m) %in% colnames(post$top_loci), info = m)
  }
  expect_false("pip" %in% colnames(post$top_loci))
})

test_that("format_finemapping_output no longer copies pip_<primary_method> into an unsuffixed pip column", {
  d <- .make_univariate_data(seed = 15, effect_idx = c(20, 40))
  fits <- fit_susie_inf_then_susie(d$X, d$y)
  post <- postprocess_finemapping_fits(fits, data_x = d$X, data_y = d$y, coverage = 0.95)
  out <- format_finemapping_output(post, primary_method = "susie")
  expect_false("pip" %in% colnames(out$top_loci))
  expect_true("pip_susie" %in% colnames(out$top_loci))
})

test_that(".translate_legacy_top_loci_cs_columns renames pip_susie -> pip for legacy callers", {
  new_format <- data.frame(
    variant_id = c("v1", "v2"),
    pip_susie = c(0.9, 0.1),
    CS_95_susie = c(1, 0),
    pip_susie_inf = c(0.8, 0.2),
    CS_95_susie_inf = c(1, 0),
    stringsAsFactors = FALSE
  )
  out <- pecotmr:::.translate_legacy_top_loci_cs_columns(new_format)
  expect_true("pip" %in% colnames(out))
  expect_false("pip_susie" %in% colnames(out))
  # The susie_inf and CS columns are untouched
  expect_true("pip_susie_inf" %in% colnames(out))
  expect_true("CS_95_susie" %in% colnames(out))
  expect_true("CS_95_susie_inf" %in% colnames(out))
})

test_that(".translate_legacy_top_loci_cs_columns leaves existing pip column alone", {
  legacy <- data.frame(
    variant_id = c("v1", "v2"),
    pip = c(0.9, 0.1),
    cs_coverage_0.95 = c(1, 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- pecotmr:::.translate_legacy_top_loci_cs_columns(legacy)
  expect_true("pip" %in% colnames(out))
  expect_true("CS_95_susie" %in% colnames(out))   # legacy cs_coverage rename
  expect_false("cs_coverage_0.95" %in% colnames(out))
  expect_false("pip_susie" %in% colnames(out))     # no double-conversion
})


test_that("top_loci has symmetric CS_<cov>_<method> columns even for methods with no surviving CSes", {
  d <- .make_univariate_data(seed = 16, effect_idx = c(25))
  fits <- fit_susie_inf_then_susie(d$X, d$y)
  # Force one method (susie_inf) to have no surviving effects: zero out V so
  # select_effects drops everything.
  fits$susie_inf$V <- rep(0, length(fits$susie_inf$V))
  post <- postprocess_finemapping_fits(fits, data_x = d$X, data_y = d$y,
                                       coverage = 0.95, secondary_coverage = c(0.7, 0.5))
  cols <- colnames(post$top_loci)
  for (m in c("susie", "susie_inf")) {
    for (cov_lbl in c("CS_95_", "CS_70_", "CS_50_")) {
      expect_true(paste0(cov_lbl, m) %in% cols,
                  info = paste("missing", paste0(cov_lbl, m)))
    }
  }
  # susie_inf had no surviving effects -> all CS_*_susie_inf are 0
  for (col in grep("^CS_.*_susie_inf$", cols, value = TRUE)) {
    expect_true(all(post$top_loci[[col]] == 0L), info = col)
  }
})

# ============================================================================
# Unified top-loci annotated long + build_top_loci_export
# ============================================================================

# Helper: build a minimal fake fit + cs_tables that build_top_loci_long accepts.
.make_fake_inputs_for_annotated_long <- function(variant_ids,
                                                 cs_membership,
                                                 cs_purity_value = 0.85,
                                                 n_samples = 100,
                                                 n_variants = NULL,
                                                 pip = NULL) {
  p <- length(variant_ids)
  if (is.null(n_variants)) n_variants <- p
  if (is.null(pip)) pip <- seq(0.6, 0.9, length.out = p)
  alpha <- matrix(0, nrow = 1, ncol = p,
                  dimnames = list(NULL, variant_ids))
  alpha[1, ] <- pip / sum(pip)
  mu <- matrix(0.5, nrow = 1, ncol = p,
               dimnames = list(NULL, variant_ids))
  mu2 <- mu^2 + 0.1
  fit <- list(pip = setNames(pip, variant_ids), alpha = alpha,
              mu = mu, mu2 = mu2)
  purity_df <- data.frame(min.abs.corr = cs_purity_value,
                          mean.abs.corr = cs_purity_value,
                          median.abs.corr = cs_purity_value)
  cs_table <- list(
    sets = list(cs = list(L1 = cs_membership),
                cs_index = 1L,
                requested_coverage = 0.95,
                purity = purity_df),
    cs_corr = list(matrix(c(1, cs_purity_value, cs_purity_value, 1),
                          nrow = 2)),
    pip = fit$pip
  )
  cs_tables <- list(cs_table)
  attr(cs_tables, "coverage") <- 0.95
  X <- matrix(0, nrow = n_samples, ncol = n_variants)
  rownames(X) <- paste0("sample", seq_len(n_samples))
  if (n_variants == length(variant_ids)) colnames(X) <- variant_ids
  Y <- matrix(0, nrow = n_samples, ncol = 1,
              dimnames = list(paste0("sample", seq_len(n_samples)),
                              "ENSG00000179403"))
  list(fit = fit, cs_tables = cs_tables, variant_names = variant_ids,
       data_x = X, data_y = Y)
}

test_that("build_top_loci_long emits annotated columns from fit + cs_tables", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .make_fake_inputs_for_annotated_long(variant_ids,
                                              cs_membership = c(1L, 2L))
  long <- pecotmr:::build_top_loci_long(
    fit = inp$fit, cs_tables = inp$cs_tables,
    variant_names = inp$variant_names,
    sumstats = list(betahat = c(0.2, -0.1),
                    sebetahat = c(0.05, 0.04),
                    z = c(4, -2.5)),
    maf = c(0.10, 0.25), method = "susie", signal_cutoff = 0.05
  )
  expect_true(all(c("conditional_effect", "conditional_effect_se",
                    "cs_purity") %in% names(long)))
  expect_equal(length(unique(long$cs_purity)), 1L)
  expect_equal(round(unique(long$cs_purity), 3), 0.85)
  expect_true(all(is.finite(long$conditional_effect)))
  expect_true(all(long$conditional_effect_se >= 0))
})

test_that("build_top_loci_long writes per-fit constants from data_x/data_y/other_quantities", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .make_fake_inputs_for_annotated_long(variant_ids,
                                              cs_membership = c(1L, 2L),
                                              n_samples = 419,
                                              n_variants = 11332)
  other_q <- list(
    dropped_samples = list(X = character(), y = character()),
    region          = "chr10:10823338-14348298",
    condition_id    = "Ast_DeJager_eQTL"
  )
  long <- pecotmr:::build_top_loci_long(
    fit = inp$fit, cs_tables = inp$cs_tables,
    variant_names = inp$variant_names,
    sumstats = list(betahat = c(0.2, -0.1),
                    sebetahat = c(0.05, 0.04),
                    z = c(4, -2.5)),
    maf = c(0.10, 0.25), method = "susie", signal_cutoff = 0.05,
    data_x = inp$data_x, data_y = inp$data_y, other_quantities = other_q
  )
  expect_equal(unique(long$n), 419L)
  expect_equal(unique(long$variant_number), 11332L)
  expect_equal(unique(long$gene_id), "ENSG00000179403")
  expect_equal(unique(long$region), "chr10:10823338-14348298")
  expect_equal(unique(long$event_ID), "Ast_DeJager_eQTL_ENSG00000179403")
})

test_that("build_top_loci_long omits region/event_ID columns when caller does not supply them", {
  variant_ids <- c("chr1:100:A:G")
  inp <- .make_fake_inputs_for_annotated_long(variant_ids,
                                              cs_membership = 1L)
  long <- pecotmr:::build_top_loci_long(
    fit = inp$fit, cs_tables = inp$cs_tables,
    variant_names = inp$variant_names,
    sumstats = list(betahat = 0.2, sebetahat = 0.05, z = 4),
    maf = 0.10, method = "susie", signal_cutoff = 0.05,
    data_x = inp$data_x, data_y = inp$data_y, other_quantities = NULL
  )
  expect_false("region" %in% names(long))
  expect_false("event_ID" %in% names(long))
  expect_true("n" %in% names(long))
  expect_true("variant_number" %in% names(long))
  expect_true("gene_id" %in% names(long))
})

test_that("build_top_loci_long stays backward-compatible when data_x/data_y/other_quantities are omitted", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .make_fake_inputs_for_annotated_long(variant_ids,
                                              cs_membership = c(1L, 2L))
  long <- pecotmr:::build_top_loci_long(
    fit = inp$fit, cs_tables = inp$cs_tables,
    variant_names = inp$variant_names,
    sumstats = list(betahat = c(0.2, -0.1),
                    sebetahat = c(0.05, 0.04),
                    z = c(4, -2.5)),
    maf = c(0.10, 0.25), method = "susie", signal_cutoff = 0.05
  )
  # n / variant_number / gene_id columns exist but are NA because no
  # data_x / data_y were supplied; region / event_ID are omitted entirely.
  expect_true(all(is.na(long$n)))
  expect_true(all(is.na(long$variant_number)))
  expect_true(all(is.na(long$gene_id)))
  expect_false("region" %in% names(long))
  expect_false("event_ID" %in% names(long))
})

test_that("build_top_loci_export errors when required columns are missing", {
  bad <- data.frame(variant_id = "chr1:100:A:G", pip = 0.9,
                    coverage = 0.95, cs = 1L,
                    stringsAsFactors = FALSE)
  expect_error(build_top_loci_export(bad),
               "missing required columns")
})

test_that("build_top_loci_export returns the fixed 21-column schema on an empty long", {
  empty <- pecotmr:::.empty_top_loci_long()
  out <- build_top_loci_export(empty)
  expected <- c("#chr", "start", "end", "a1", "a2", "variant_ID", "gene_ID",
                "event_ID", "cs_coverage_0.95", "cs_coverage_0.7",
                "cs_coverage_0.5", "cs_purity", "PIP", "conditional_effect",
                "conditional_effect_se", "analysis_region",
                "analysis_variants_number", "beta", "se", "n", "maf")
  expect_equal(names(out), expected)
  expect_equal(nrow(out), 0L)
})

test_that("build_top_loci_export projects an annotated long into the fixed compact schema", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .make_fake_inputs_for_annotated_long(variant_ids,
                                              cs_membership = c(1L, 2L),
                                              n_samples = 419,
                                              n_variants = 11332)
  other_q <- list(region = "chr10:10823338-14348298",
                  condition_id = "Ast_DeJager_eQTL")
  long <- pecotmr:::build_top_loci_long(
    fit = inp$fit, cs_tables = inp$cs_tables,
    variant_names = inp$variant_names,
    sumstats = list(betahat = c(0.2, -0.1),
                    sebetahat = c(0.05, 0.04),
                    z = c(4, -2.5)),
    maf = c(0.10, 0.25), method = "susie", signal_cutoff = 0.05,
    data_x = inp$data_x, data_y = inp$data_y, other_quantities = other_q
  )
  out <- build_top_loci_export(long)
  expected_cols <- c("#chr", "start", "end", "a1", "a2", "variant_ID",
                     "gene_ID", "event_ID", "cs_coverage_0.95",
                     "cs_coverage_0.7", "cs_coverage_0.5", "cs_purity",
                     "PIP", "conditional_effect", "conditional_effect_se",
                     "analysis_region", "analysis_variants_number",
                     "beta", "se", "n", "maf")
  expect_equal(names(out), expected_cols)
  expect_equal(out$gene_ID[[1]], "ENSG00000179403")
  expect_equal(out$event_ID[[1]], "Ast_DeJager_eQTL_ENSG00000179403")
  expect_equal(out$analysis_region[[1]], "chr10:10823338-14348298")
  expect_equal(out$analysis_variants_number[[1]], 11332L)
  expect_equal(out$n[[1]], 419L)
})
