# krls -- Kernel Regularized Least Squares.
#
# File layout:
#   1. krls()              -- S3 generic dispatch.
#   2. krls.formula()      -- formula -> matrix dispatch (mirrors ares.formula).
#   3. krls.default()      -- main fitter, mirrors KRLS::krls (Hainmueller and
#                             Hazlett 2014) numerically.  Returns an S3
#                             object of class c("krls_rr", "krls").
#   4. predict.krls_rr()   -- predictions on new data, optional fit SEs.
#                             Handles formula-side terms / factor expansion.
#   5. print.krls_rr()     -- compact one-screen summary.
#   6. summary.krls_rr()   -- richer summary with quartiles of marginal
#                             effects (mirrors KRLS::summary.krls layout).
#   7. .krls_lambdasearch  -- golden-section LOO lambda selector.
#   8. .krls_fd_binary     -- replace binary-column derivatives with
#                             finite-difference contrasts (mirrors
#                             KRLS::fdskrls).
#
# Conventions:
#   - Functions prefixed with `.` are internal helpers (not exported).
#   - The C++ engine is reached via Rcpp exports: krls_kernel_cpp,
#     krls_kernel_pred_cpp, krls_eig_cpp, krls_solve_cpp,
#     krls_loo_loss_cpp, krls_deriv_cpp, krls_avg_deriv_var_cpp,
#     krls_vsq_cpp (see src/krls.cpp).
#   - Output field names follow KRLS::krls (`coeffs`, `Looe`, `fitted`,
#     `sigma`, `lambda`, `R2`, `derivatives`, `avgderivatives`,
#     `var.avgderivatives`, `vcov.c`, `vcov.fitted`, `binaryindicator`)
#     so downstream code that consumed KRLS output works unchanged.
#
# Naming note: KRLS exports its own krls() with class "krls".  Loading
# both packages masks one or the other depending on attach order.  Use
# `roadrunner::krls()` or `KRLS::krls()` to disambiguate.

#' Kernel Regularized Least Squares
#'
#' Fits a Kernel Regularized Least Squares model with Gaussian kernel,
#' selecting the ridge penalty by leave-one-out cross-validation via the
#' closed-form identity of Hainmueller and Hazlett (2014).  Marginal
#' effects and their variances are computed by default and returned on
#' the original `(X, y)` scale.
#'
#' Three call styles are supported (mirrors `ares()`):
#'
#'   * `krls(X, y, ...)`  -- matrix / numeric interface (back-compatible).
#'   * `krls(y ~ x1 + x2 + ..., data = df, ...)` -- formula interface
#'     with factor expansion and derived terms (`I(x^2)`, `poly(x, 2)`,
#'     etc).
#'   * `krls(df, y, ...)` where `df` is a data frame with numeric +
#'     factor / character columns -- expands categorical columns via
#'     `model.matrix(~ ., df)` (treatment contrasts, intercept dropped).
#'
#' The numerical pipeline mirrors `KRLS::krls()` exactly:
#'
#' 1. `X` and `y` are standardised (column-centred, unit sd).
#' 2. The Gaussian kernel `K_ij = exp(-||x_i - x_j||^2 / sigma)` is built
#'    in parallel C++.
#' 3. `K` is eigendecomposed.  All subsequent solves use the
#'    eigen-basis closed forms, never inverting `K + lambda I` directly.
#' 4. `lambda` is selected by golden-section search on the closed-form
#'    LOO error sum, with the same `(L, U, tol)` bracket as `KRLS::krls`.
#' 5. Marginal effects use the closed-form identity
#'    `dy/dx_k = -(2/sigma) * (X_k * (K c) - K diag(c) X)_k`, computed
#'    without forming the `n x n` distance matrix.
#' 6. Output is unstandardised back to the original `(X, y)` scale.
#'
#' At a fixed `(sigma, lambda)`, fits agree with `KRLS::krls()` to
#' floating-point precision (typically `< 1e-12` on coefficients,
#' fitted values, and marginal effects for `n <= 1000`).  Wall-clock
#' time is roughly `6-10x` faster than `KRLS::krls()` at `n >= 500`
#' when marginal effects and variance estimates are requested.
#'
#' Memory scales as `O(n^2)`: the kernel and its squared eigenvector
#' matrix are both stored.  Expect about `0.4 * n^2 / 1e6` MB of
#' peak working memory (e.g. ~400MB at `n = 1000`, ~10GB at `n = 5000`).
#'
#' @param X A numeric matrix or data frame of predictors (`n x p`).
#'   Constant columns are rejected. Factor / character columns in a
#'   data frame are expanded via `model.matrix(~ ., x)` (treatment
#'   contrasts, intercept dropped). Missing values are handled via
#'   `na.action` (see below).
#' @param y A numeric response vector or single-column matrix.  Constant
#'   `y` is rejected.
#' @param data Used only by the formula method. A data frame containing
#'   the variables referenced by the formula.
#' @param subset Used only by the formula method. An optional integer
#'   or logical vector restricting rows of `data` used for the fit.
#' @param sigma Gaussian-kernel bandwidth.  Default `NULL`, which sets
#'   sigma via the geomean_p formula: `sqrt(median(d2) * p)` where `d2`
#'   are pairwise squared Euclidean distances on the standardised
#'   predictors and `p = ncol(X)`.  Must be a positive scalar if supplied.
#' @param lambda Optional ridge penalty.  If `NULL` (default), selected
#'   by golden-section search on the LOO error.
#' @param derivative Logical.  If `TRUE` (default), compute pointwise
#'   marginal effects and their average per variable.  Requires `vcov`.
#' @param binary Logical.  If `TRUE` (default), columns of `X` with
#'   exactly two unique values are treated as binary and their marginal
#'   effects are replaced by predicted-Y first differences (matches
#'   `KRLS::fdskrls`).
#' @param vcov Logical.  If `TRUE` (default), compute the coefficient
#'   covariance and the variance of average marginal effects.
#' @param weights Optional vector of observation weights (length `n`,
#'   strictly positive). Internally normalised to mean 1. Implements
#'   weighted KRLS via a `D K D` transform where `D = diag(sqrt(w))`.
#'   `weights = rep(1, n)` is byte-identical to the unweighted path.
#' @param L,U Optional lower / upper bracket for the lambda search.  If
#'   `NULL`, defaults follow `KRLS::krls()`.
#' @param tol Tolerance for the lambda golden section.  Default `1e-6`
#'   (fixed, independent of `n`; was `1e-3 * n` in earlier versions).
#'   A fixed small tolerance gives 6-digit lambda precision regardless
#'   of sample size.
#' @param eigtrunc Optional eigenvalue truncation cutoff in `(0, 1]`.
#'   When set, eigenvalues below `eigtrunc * max(d)` are dropped from
#'   the solve.  `NULL` (default) keeps all eigenvalues.
#' @param lambda.method Lambda-selection rule. `"loo"` (default) uses
#'   the closed-form leave-one-out golden-section search; `"cv"`
#'   uses K-fold CV over a grid (`nfold > 0` required). `"gcv"` uses
#'   the closed-form generalised cross-validation criterion (Craven &
#'   Wahba 1979) using the same eigendecomposition; recommended when
#'   `n` is large or LOO behaves unstably.
#' @param lambda.grid Optional numeric vector of lambda candidates for
#'   `lambda.method = "cv"`. `NULL` (default) auto-generates a
#'   log-spaced grid in `[L, U]`.
#' @param nfold Number of CV folds for `lambda.method = "cv"` or for
#'   `autotune`. `0` (default) disables CV.
#' @param ncross Number of CV repetitions (each builds a fresh fold
#'   partition). Default `NULL`, which resolves to `2` for the autotune
#'   sigma search (repeated-CV stabilisation) and `1` for
#'   `lambda.method = 'cv'`. Explicit integer values are honoured
#'   in both paths.
#' @param stratify If `TRUE` (default), CV folds are quantile-
#'   stratified on `y`.
#' @param seed.cv Optional integer seed for the CV fold partition.
#' @param cv.1se If `TRUE`, applies the one-standard-error rule when
#'   picking lambda under CV (smallest model within 1 SE of the
#'   minimum mean CV-MSE). Default `FALSE`.
#' @param autotune If `TRUE`, runs an inner repeated-CV grid search
#'   over `sigma` (default: 9-point multiplicative grid centred on the
#'   median-heuristic sigma anchor, from `anchor * 0.125` to
#'   `anchor * 32`) and refits the winner on the full data using the
#'   one-standard-error rule (selects the largest sigma within 1 SE of
#'   the minimum CV-MSE). Default `FALSE`.
#' @param autotune.grid Optional numeric vector of `sigma` candidates
#'   for autotune. `NULL` (default) uses the 9-point anchor-centred
#'   grid `sigma_anchor * c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)`
#'   where `sigma_anchor` is the geomean_p anchor
#'   (`sqrt(median(d2) * p)` on the standardised predictors).
#' @param autotune.nthreads Integer. Number of worker threads used to
#'   evaluate the autotune `sigma` grid in parallel. `NULL` (default)
#'   reads from `getOption("roadrunner.krls.autotune.nthreads")`, and
#'   falls back to `1L` if unset. Capped at `length(autotune.grid)`.
#'   Only used when `autotune = TRUE`.
#' @param varmod Residual variance model used to construct prediction
#'   intervals via `predict(..., interval = "pint")`. `"none"`
#'   (default) disables PIs; `"const"` uses a homoscedastic
#'   `sigma_hat` estimated from the training residuals.
#' @param n.boot Number of bootstrap replicates for bagging. `0`
#'   (default) disables bagging. When `n.boot > 0`, prediction
#'   averages across `n.boot` replicate fits.
#' @param na.action How to handle missing values in `X`. `"impute"`
#'   (default) replaces NAs with the column median (stored on the fit
#'   and reapplied at `predict()` time). `"omit"` drops rows with any
#'   NA. Missing `y` is always an error.
#' @param trace Integer. `0` is silent, `> 0` enables progress
#'   diagnostics (currently: prints chosen lambda when `> 1`,
#'   golden-section progress when `> 2`). Replaces `print.level`.
#' @param nthreads Integer. Number of threads to use for the C++ kernel
#'   build / decomposition. `0` (default) means use
#'   `RcppParallel::defaultNumThreads()`.
#' @param print.level Deprecated alias for `trace`. Prefer `trace`.
#' @param approx Either `"exact"` (default) or `"nystrom"`. The exact
#'   path runs the original O(n^3) eigendecomposition on the full kernel.
#'   `"nystrom"` replaces it with a low-rank Nystrom approximation
#'   anchored at `nystrom_m` landmarks, giving O(m^3) + O(n m^2) cost.
#'   At default `nystrom_m = ceiling(sqrt(n) * 3)` test RMSE is typically
#'   within +10% of exact on smooth DGPs.
#' @param nystrom_m Integer. Number of landmarks when `approx = "nystrom"`.
#'   Default `NULL` resolves to `ceiling(sqrt(n) * 3)` (e.g. n=2000 -> m=135,
#'   n=5000 -> m=213). Tuned empirically to keep test-set RMSE within +10%
#'   of exact on smooth DGPs.
#'   Ignored when `landmarks` is supplied explicitly.
#' @param landmarks NULL (default; auto-draw via `landmark_method`), an
#'   integer vector of row indices into X, or an m x d numeric matrix of
#'   landmark coordinates in the original X scale (auto-standardized
#'   internally using the training centers/scales).
#' @param landmark_method Either `"random"` (default; uniform subsample of
#'   training rows) or `"kmeans"` (Hartigan-Wong centers on standardized X;
#'   not bit-stable across R versions). Ignored when `landmarks` is supplied.
#' @param landmark_seed Optional integer seed for the local landmark draw.
#'   When supplied, uses `.with_seed()` to avoid disturbing the caller's
#'   global RNG state. Required for byte-identical fits across consecutive
#'   `krls(approx="nystrom")` calls.
#' @param nystrom_eps Numeric. Relative ridge floor applied to the
#'   landmark-kernel eigenvalues: `D_reg = max(D, nystrom_eps * max(D))`.
#'   Defaults to `1e-9`. Stabilizes the m x m eigendecomposition when
#'   landmarks are near-collinear; the fit object's
#'   `nystrom_diagnostics$floored_count` reports how many eigenvalues hit
#'   this floor (useful for tuning m).
#' @param ... Currently unused (caught for forward compatibility).
#'
#' @details
#' **Scale-aware sigma default (geomean_p anchor, v0.0.0.9041)**
#'
#' When `sigma = NULL`, roadrunner sets the Gaussian-kernel bandwidth using the
#' *geomean_p* formula: let `d2_ij = ||Xs_i - Xs_j||^2` be the pairwise squared
#' Euclidean distances on the standardised predictor matrix `Xs`; then
#' `sigma = sqrt(median({d2_ij : i < j}) * ncol(Xs))`.  This is the geometric
#' mean of the raw median heuristic and `p = ncol(X)`, and empirically lands in
#' the sweet spot between the two extremes.  At `n=500, p=10` on standardised
#' N(0,1) data it yields `sigma ~ 13.5`, versus `~ 18-19` for the raw median
#' and `10` for `KRLS::krls()`'s fixed default.
#'
#' **Why geomean_p**: a 15-DGP head-to-head vs `KRLS::krls()` (REQ-20260518-003,
#' iter-2) with the geomean_p anchor shows roadrunner **wins 12/15 DGPs, loses
#' 0/15, ties 3/15**.  The raw median anchor (v0.0.0.9040) won only 4/15 DGPs
#' because it over-smoothed locally nonlinear signals (exp-decay, interaction,
#' tanh-interaction, poly2) by selecting sigma ~2x too wide.  The geomean_p
#' anchor preserves scale-awareness (adapts to the actual data distribution,
#' unlike the hard-coded `sigma = p`) while avoiding over-smoothing.
#'
#' For `n > 500` the pairwise distance matrix is `O(n^2)`; to keep the default
#' cheap, a 500-row subsample is drawn with a fixed seed (`set.seed(2718)`) so
#' the result is deterministic within a session.
#'
#' **Autotune-equals-default equivalence at well-fit settings**
#'
#' When autotune is enabled at settings where the default sigma is already near-
#' optimal (e.g. additive and interaction DGPs at `n=500, p=10`), the autotune
#' grid is centred on `sigma_anchor` and the CV argmin coincides with the grid
#' centre.  The 1-SE rule then selects `sigma_anchor`, producing a fit identical
#' to the non-autotuned default.  This is the expected behaviour — it confirms
#' that the default sigma is well-chosen for these DGPs — not a failure of
#' autotune.  Autotune yields improvements when the optimal sigma departs from
#' the anchor (e.g. sparse / linear DGPs where a much wider kernel is better).
#'
#' Since v0.0.0.9042 the autotune inner loop is parallelised over sigma
#' candidates using `RcppParallel` (TBB). The pairwise squared-distance
#' matrix is computed once per CV fold and reused across every sigma in
#' the grid (`exp(-D / sigma)` is then a cheap elementwise op).
#' Determinism is preserved: each worker writes to a unique slot in the
#' output vectors, so the result is byte-identical regardless of
#' `autotune.nthreads`. Pass `autotune.nthreads = 1` to force strictly
#' sequential execution.
#'
#' @note
#' At a fixed `(sigma, lambda)` fits remain byte-identical to all earlier
#' versions of roadrunner KRLS.  The four defaults changed in v0.0.0.9040
#' (`sigma`, `tol`, autotune nfold/ncross, autotune grid) only affect results
#' when those arguments are left at `NULL` / at their default.  Restore any
#' old default explicitly: `sigma = ncol(X)`, `tol = 1e-3 * nrow(X)`,
#' `nfold = 5`, `ncross = 1`, `autotune.grid = ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)`.
#'
#' The sigma anchor formula was refined in v0.0.0.9041 from the raw median
#' (`median(d2)`) to the geomean_p formula (`sqrt(median(d2) * p)`).  To
#' restore the v0.0.0.9040 median anchor, pass
#' `sigma = stats::median(as.numeric(stats::dist(scale(X)))^2)`.
#'
#' Phase 1 speedup (v0.0.0.9042): parallel autotune assumes single-threaded
#' BLAS for the small kernel operations inside the parallel region. If you
#' have set a high process-global BLAS thread count (e.g. via
#' `RhpcBLASctl::blas_set_num_threads(8)`), consider resetting to 1 around
#' `krls(..., autotune = TRUE)` calls to avoid oversubscription. The
#' BLAS-heavy distance computation runs OUTSIDE the parallel region and
#' benefits from multi-threaded BLAS.
#'
#' `approx = "nystrom"` is currently incompatible with observation
#' `weights` and with user-supplied lambda-search bracket `L`/`U`/`tol`.
#' Both error with a clear message at fit time. `vcov = TRUE` is
#' supported in the single-fit Nystrom path but the resulting `vcov`
#' is for the m-length dual coefficients, not the n-length kernel
#' coefficients. Phase 3 may extend these.
#'
#' Since v0.0.0.9043 `krls()` supports an opt-in Nystrom low-rank
#' approximation via `approx = "nystrom"`. Replacing the full n x n
#' eigendecomposition with an m x m one (m = `nystrom_m`, default
#' `ceiling(sqrt(n) * 3)`) gives ~5x speedup at n=2000 and ~10x at
#' n=5000 with RMSE typically within +10% of the exact fit on smooth
#' DGPs. Landmarks default to a uniform random subsample of training
#' rows; pass `landmark_method = "kmeans"` for centroidal landmarks.
#' Fits are byte-identical at fixed `landmark_seed` and `seed.cv`.
#'
#' @return An object of S3 class `c("krls_rr", "krls")` with components
#'   mirroring `KRLS::krls()`: `K`, `coeffs`, `Looe`, `fitted`, `X`,
#'   `y`, `sigma`, `lambda`, `R2`, `derivatives`, `avgderivatives`,
#'   `var.avgderivatives`, `vcov.c`, `vcov.fitted`, `binaryindicator`.
#'   Formula-method fits additionally carry `call`, `terms`, `xlevels`,
#'   `factor_info`, `na.action`, and `na.medians` for downstream
#'   `predict()`, `update()`, and `model.matrix()` support.
#'
#'   `Looe` follows the `KRLS::krls()` scale convention: it is the sum
#'   of squared leave-one-out residuals on the *standardised* `y` scale
#'   multiplied by `sd(y)`, so its units are `[y^2 / sd_y] = [y]`.  This
#'   is preserved for downstream compatibility with code that consumed
#'   `KRLS::krls()` output; it is **not** the LOO MSE in raw-`y`
#'   squared units.
#'
#' @references Hainmueller, J. and C. Hazlett (2014).  "Kernel
#'   Regularized Least Squares: Reducing Misspecification Bias with a
#'   Flexible and Interpretable Machine Learning Approach."  *Political
#'   Analysis* 22(2):143--168.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(rnorm(n * 3), n, 3)
#' colnames(X) <- c("age", "income", "score")
#' y <- sin(X[, 1]) + 0.5 * X[, 2]^2 - 0.3 * X[, 3] +
#'   rnorm(n, sd = 0.2)
#'
#' fit <- krls(X, y)
#' fit
#' fit$avgderivatives           # average marginal effect per variable
#' summary(fit)
#'
#' ## Formula interface with mixed numeric + factor predictors.
#' df <- data.frame(age = X[, 1], income = X[, 2],
#'                  region = factor(sample(c("N", "S"), n, replace = TRUE)),
#'                  y = y)
#' fit_f <- krls(y ~ age + income + region, data = df)
#'
#' ## Predictions on new data with pointwise SEs.
#' Xnew <- matrix(rnorm(20 * 3), 20, 3)
#' colnames(Xnew) <- colnames(X)
#' pr <- predict(fit, Xnew, se.fit = TRUE)
#' head(pr$fit)
#' head(pr$se.fit)
#'
#' @export
krls <- function(X, ...) UseMethod("krls")

#' @rdname krls
#' @export
krls.formula <- function(X, data = NULL, subset = NULL, ..., y = NULL) {
  formula <- X
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()

  # Reject offset() terms: KRLS has no GLM offset concept; absorbing into
  # the design silently would be wrong.
  pre_tt <- stats::terms(formula, data = data)
  if (length(attr(pre_tt, "offset"))) {
    stop("krls: offset() terms are not supported in the formula. ",
         "Subtract the offset from y before calling krls() if needed.",
         call. = FALSE)
  }

  # Build model frame with na.pass so krls.default applies its own
  # na.action. Subset is honoured via the lm()-style match.call() trick
  # so NSE resolves it against the caller's environment.
  mfcall <- match.call(expand.dots = FALSE)
  mfcall$y <- NULL
  mfcall[[1L]] <- quote(stats::model.frame)
  names(mfcall)[names(mfcall) == "X"] <- "formula"
  mfcall$`...` <- NULL
  mfcall$na.action <- quote(stats::na.pass)
  mf <- eval(mfcall, parent.frame())
  yv <- stats::model.response(mf)
  if (is.null(yv))
    stop("krls: response variable is missing from formula/data.")
  mm <- stats::model.matrix(formula, mf)
  has_int <- "(Intercept)" %in% colnames(mm)
  if (has_int) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  out <- krls.default(X = mm, y = as.numeric(yv), ...)
  out$call <- cl
  tt <- stats::terms(formula, data = data)
  out$terms <- tt
  out$xlevels <- stats::.getXlevels(tt, mf)
  out
}

#' @rdname krls
#' @export
krls.default <- function(X, y,
                         sigma = NULL, lambda = NULL,
                         derivative = TRUE, binary = TRUE, vcov = TRUE,
                         weights = NULL,
                         L = NULL, U = NULL, tol = NULL, eigtrunc = NULL,
                         lambda.method = c("loo", "gcv", "cv"),
                         lambda.grid = NULL,
                         nfold = 0L, ncross = NULL, stratify = TRUE,
                         seed.cv = NULL, cv.1se = FALSE,
                         autotune = FALSE, autotune.grid = NULL,
                         autotune.nthreads = NULL,
                         varmod = c("none", "const"),
                         n.boot = 0L,
                         na.action = c("impute", "omit"),
                         approx = c("exact", "nystrom"),
                         nystrom_m = NULL,
                         landmarks = NULL,
                         landmark_method = c("random", "kmeans"),
                         landmark_seed = NULL,
                         nystrom_eps = 1e-9,
                         trace = NULL, nthreads = 0L,
                         print.level = NULL, ...) {
  cl <- match.call()
  na.action <- match.arg(na.action)
  lambda.method <- match.arg(lambda.method)
  varmod <- match.arg(varmod)
  approx <- match.arg(approx)
  landmark_method <- match.arg(landmark_method)
  nystrom_eps <- .validate_nystrom_eps(nystrom_eps)

  if (lambda.method == "gcv" && approx == "nystrom") {
    stop("krls: lambda.method = 'gcv' is not supported with approx = 'nystrom'.")
  }

  ## --- print.level / trace harmonisation ---------------------------
  # Phase 8: `trace` is the canonical name (matches ares).  `print.level`
  # is kept as a deprecated alias for back-compat; we warn once if the
  # user supplies it AND it would shadow an explicit `trace`.
  if (!is.null(print.level)) {
    if (is.null(trace)) {
      trace <- print.level
    } else if (!identical(as.integer(trace), as.integer(print.level))) {
      warning("krls: both `trace` and `print.level` supplied; using `trace`.",
              " (`print.level` is deprecated; prefer `trace`.)",
              call. = FALSE)
    } else {
      warning("krls: `print.level` is deprecated; use `trace` instead.",
              call. = FALSE)
    }
  }
  if (is.null(trace)) trace <- 0L
  trace <- as.integer(trace)

  ## --- argument validation -----------------------------------------
  if (is.null(X) || is.null(y)) stop("X and y are required")
  if (is.factor(y) || is.character(y)) {
    stop("y must be numeric (got ",
         if (is.factor(y)) "factor" else "character", ")")
  }
  if (is.factor(X) || is.character(X)) {
    stop("X must be numeric (got ",
         if (is.factor(X)) "factor" else "character", ")")
  }

  ## --- data.frame factor / character expansion ---------------------
  factor_info <- NULL
  pre_na_medians <- NULL
  if (is.data.frame(X)) {
    is_cat <- vapply(X, function(z) is.factor(z) || is.character(z),
                     logical(1L))
    is_num <- vapply(X, is.numeric, logical(1L))

    # Numeric-column NA handling BEFORE model.matrix (which would drop
    # rows under its default na.action). Mirrors ares.default.
    if (any(is_num)) {
      num_names <- names(X)[is_num]
      df_na_mask <- vapply(num_names,
                           function(jn) anyNA(X[[jn]]), logical(1L))
      if (any(df_na_mask)) {
        rows_aff <- sum(rowSums(is.na(
          as.matrix(X[, num_names, drop = FALSE]))) > 0L)
        cols_aff <- sum(df_na_mask)
        if (na.action == "impute") {
          pre_na_medians <- vapply(
            num_names,
            function(jn) stats::median(X[[jn]], na.rm = TRUE),
            numeric(1L))
          names(pre_na_medians) <- num_names
          if (any(!is.finite(pre_na_medians)))
            stop("krls: na.action='impute' but at least one numeric column",
                 " is entirely NA; cannot compute median.",
                 " Drop the column or pass na.action='omit'.")
          for (jn in num_names) {
            na_j <- is.na(X[[jn]])
            if (any(na_j)) X[[jn]][na_j] <- pre_na_medians[[jn]]
          }
          warning("krls: median-imputed ", sum(df_na_mask),
                  " numeric column(s) across ", rows_aff,
                  " row(s). Column medians are stored on the fit.",
                  call. = FALSE)
        } else {
          keep <- stats::complete.cases(X[, num_names, drop = FALSE])
          dropped <- sum(!keep)
          warning("krls: dropped ", dropped,
                  " incomplete row(s) across ", cols_aff,
                  " numeric column(s). Pass na.action='impute' to",
                  " median-impute instead.", call. = FALSE)
          X <- X[keep, , drop = FALSE]
          y <- y[keep]
          if (!is.null(weights)) weights <- weights[keep]
          if (length(y) < 3L)
            stop("krls: na.action='omit' left fewer than 3 rows; aborting.")
        }
      }
    }

    if (any(is_cat)) {
      for (j in which(is_cat)) {
        if (is.character(X[[j]])) X[[j]] <- factor(X[[j]])
      }
      xlevels <- lapply(X[is_cat], levels)
      xmm <- stats::model.matrix(~ ., data = X)
      if ("(Intercept)" %in% colnames(xmm))
        xmm <- xmm[, colnames(xmm) != "(Intercept)", drop = FALSE]
      factor_info <- list(
        xlevels       = xlevels,
        orig_names    = names(X),
        is_cat        = is_cat,
        is_num        = is_num,
        expanded_names = colnames(xmm),
        num_medians   = pre_na_medians
      )
      X <- xmm
    } else {
      nm <- names(X)
      X <- as.matrix(X)
      storage.mode(X) <- "double"
      colnames(X) <- nm
    }
  }

  X <- as.matrix(X)
  y <- as.matrix(y)
  if (ncol(y) != 1L) {
    stop("y must be a vector or single-column matrix (got ", ncol(y),
         " columns)")
  }
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  if (!is.numeric(X)) stop("X must be numeric")
  if (!is.numeric(y)) stop("y must be numeric")
  if (anyNA(y))      stop("krls: y contains missing data")

  ## --- na.action on numeric matrix ---------------------------------
  na_medians <- NULL
  if (any(is.nan(X)) || any(is.infinite(X)))
    stop("krls: X contains NaN or +/-Inf values; aborting fit.",
         " (NA is handled via `na.action`; NaN / Inf are not.)")
  na_x <- is.na(X)
  if (any(na_x)) {
    rows_aff <- sum(rowSums(na_x) > 0L)
    cols_aff <- sum(colSums(na_x) > 0L)
    if (na.action == "impute") {
      na_medians <- vapply(seq_len(ncol(X)),
                           function(j) stats::median(X[, j], na.rm = TRUE),
                           numeric(1L))
      names(na_medians) <- colnames(X)
      if (any(!is.finite(na_medians)))
        stop("krls: na.action='impute' but at least one X column is",
             " entirely NA; cannot compute median.")
      for (j in seq_len(ncol(X))) {
        na_j <- na_x[, j]
        if (any(na_j)) X[na_j, j] <- na_medians[j]
      }
      warning("krls: median-imputed ", sum(na_x), " NA value(s) across ",
              rows_aff, " row(s) and ", cols_aff, " column(s).",
              " Column medians stored on fit.", call. = FALSE)
    } else {
      keep <- rowSums(na_x) == 0L
      dropped <- sum(!keep)
      warning("krls: dropped ", dropped, " incomplete row(s) across ",
              cols_aff, " column(s).", call. = FALSE)
      X <- X[keep, , drop = FALSE]
      y <- y[keep, , drop = FALSE]
      if (!is.null(weights)) weights <- weights[keep]
      if (length(y) < 3L)
        stop("krls: na.action='omit' left fewer than 3 rows; aborting.")
    }
  }

  if (var(as.vector(y)) == 0) stop("y is a constant (does not vary)")
  stopifnot(is.logical(derivative), is.logical(vcov), is.logical(binary))
  if (derivative && !vcov) {
    stop("derivative = TRUE requires vcov = TRUE")
  }
  n <- nrow(X); d <- ncol(X)
  if (n != nrow(y)) stop("nrow(X) not equal to length of y")

  ## --- weights validation (Phase 2) -------------------------------
  w_norm <- NULL
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("krls: weights must be numeric.")
    if (length(weights) != n)
      stop("krls: length(weights) (", length(weights),
           ") must equal nrow(X) (", n, ").")
    if (any(!is.finite(weights)))
      stop("krls: weights contain NA / NaN / Inf; all weights must be finite.")
    if (any(weights <= 0))
      stop("krls: weights must be strictly positive.")
    wm <- mean(weights)
    if (wm <= 0)
      stop("krls: weights mean is zero; cannot normalise.")
    w_norm <- as.numeric(weights) / wm
  }

  if (!is.null(eigtrunc)) {
    stopifnot(is.numeric(eigtrunc), length(eigtrunc) == 1)
    if (eigtrunc < 0 || eigtrunc > 1) {
      stop("eigtrunc must be in [0, 1]")
    }
    if (eigtrunc == 0) {
      eigtrunc <- NULL
      warning("eigtrunc = 0 is equivalent to NULL; ignoring")
    }
  }
  ## sigma and autotune.grid defaults are set AFTER standardisation
  ## (they need Xs for the scale-aware sigma anchor; see Fix 1 + Fix 4).
  ## Capture whether the user supplied these args so we know whether to
  ## compute sigma_anchor.
  sigma_user_null     <- is.null(sigma)
  autotune_grid_null  <- is.null(autotune.grid)
  if (!is.null(sigma)) {
    stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  }
  if (!autotune_grid_null) {
    stopifnot(is.numeric(autotune.grid), length(autotune.grid) >= 1L,
              all(autotune.grid > 0))
    autotune.grid <- sort(unique(as.numeric(autotune.grid)))
  }

  # autotune.nthreads -- default via getOption fallback
  if (is.null(autotune.nthreads)) {
    autotune.nthreads <- getOption(
      "roadrunner.nthreads",
      max(1L, parallel::detectCores(logical = FALSE))
    )
  }
  stopifnot(is.numeric(autotune.nthreads), length(autotune.nthreads) == 1L,
            autotune.nthreads >= 1)
  autotune.nthreads <- as.integer(autotune.nthreads)
  if (length(autotune.grid) > 0L) {
    autotune.nthreads <- min(autotune.nthreads, length(autotune.grid))
  }
  if (isTRUE(RcppParallel::defaultNumThreads() == 1L)) {
    autotune.nthreads <- 1L
  }

  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(d))

  ## --- nthreads (Phase 8) -----------------------------------------
  nthreads <- as.integer(nthreads)
  if (is.na(nthreads) || nthreads < 0L) nthreads <- 0L
  nthreads_eff <- if (nthreads <= 0L)
    RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  ## --- standardise -------------------------------------------------
  X.init     <- X
  X.init.sd  <- apply(X.init, 2L, sd)
  if (any(X.init.sd == 0)) {
    stop("at least one column in X is a constant, please remove the constant(s)")
  }
  y.init      <- y
  y.init.sd   <- apply(y.init, 2L, sd)
  y.init.mean <- mean(y.init)
  Xs <- scale(X.init, center = TRUE, scale = X.init.sd)
  ys <- scale(y.init, center = y.init.mean, scale = y.init.sd)
  ## drop attributes for clean C++ handoff
  Xs <- matrix(Xs, n, d, dimnames = list(NULL, colnames(X)))
  ys <- as.numeric(ys)

  ## --- scale-aware sigma anchor (Fix 1 + Fix 4) --------------------
  ## Compute sigma_anchor (median pairwise sq. dist on Xs) whenever
  ## (a) user did not supply sigma, or (b) autotune is on and user did
  ## not supply autotune.grid. Skip if neither condition holds.
  sigma_anchor <- NULL
  if (sigma_user_null || (isTRUE(autotune) && autotune_grid_null)) {
    sigma_anchor <- .krls_sigma_anchor(Xs)
  }

  ## --- sigma default assignment (Fix 1) ----------------------------
  if (sigma_user_null) {
    sigma <- sigma_anchor
  }
  ## (non-NULL sigma already validated above)

  ## --- autotune grid default assignment (Fix 4) --------------------
  if (isTRUE(autotune)) {
    if (autotune_grid_null) {
      ## Centre the 9-point grid on sigma_anchor; 2x multiplicative spacing.
      anchor_g <- if (!is.null(sigma_anchor)) sigma_anchor else sigma
      autotune.grid <- anchor_g * c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
    }
    ## User-supplied autotune.grid already validated + sorted above.
  }

  ## --- weighted-LS transform (Phase 2) ------------------------------
  # Solve a *weighted* KRLS by applying D = diag(sqrt(w)) on both sides
  # of K. The closed-form fit on D K D with target D ys recovers the
  # weighted normal equations
  #   (K W K + lambda K) c = K W y   <=>   (K + lambda K^-1) c = W y
  # in the standardised space. Coefficients are returned in the
  # original (un-D-twisted) basis: c_final = D c_solve so that
  # K %*% c_final == K_w %*% c_solve and yhat is on the un-weighted
  # scale (n x 1 prediction on the training X). When weights = rep(1, n)
  # this collapses byte-identically to the unweighted path (sqw is the
  # vector of ones; D is the identity).
  if (!is.null(w_norm)) {
    sqw <- sqrt(w_norm)
  } else {
    sqw <- NULL
  }

  ## --- autotune over sigma (Fix 3: repeated CV + 1-SE rule) --------
  # Inner-CV grid search: for each sigma candidate, run a K-fold
  # held-out MSE scan (using LOO closed-form per-fold lambda for each
  # cell) and pick the sigma with the lowest mean held-out MSE.
  # Fix 3a: default nfold for autotune raised to 10 (from 5).
  # Fix 3b: ncross_at outer repetitions (default 2 when ncross=NULL).
  # Fix 3c: 1-SE rule selects the largest sigma within 1 SE of min MSE.
  autotune_info <- NULL
  if (isTRUE(autotune) && identical(approx, "nystrom")) {
    ## Phase 2 T5: Nystrom autotune path. Parallel sigma sweep per fold via
    ## krls_nystrom_autotune_inner_cpp; final refit on full Xs at the chosen
    ## sigma. Weights are not yet supported in the Nystrom path.
    if (!is.null(w_norm)) {
      stop("krls: weights are not supported with approx = 'nystrom'; ",
           "use approx = 'exact'", call. = FALSE)
    }
    if (!is.null(L) || !is.null(U) || !is.null(tol)) {
      stop("krls: user-supplied L/U/tol not supported with approx = 'nystrom'",
           call. = FALSE)
    }
    nfold_at <- as.integer(nfold)
    if (is.na(nfold_at) || nfold_at <= 0L) nfold_at <- 10L
    if (nfold_at > nrow(X))
      stop("krls: autotune nfold (", nfold_at,
           ") exceeds nrow(X) (", nrow(X), ").")
    ncross_at <- if (is.null(ncross)) 2L else as.integer(ncross)
    if (is.na(ncross_at) || ncross_at < 1L) ncross_at <- 1L

    if (!is.null(seed.cv)) {
      if (exists(".Random.seed", envir = globalenv(),
                 inherits = FALSE)) {
        old_seed_at <- get(".Random.seed", envir = globalenv(),
                            inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed_at,
                       envir = globalenv()), add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()),
                add = TRUE)
      }
      set.seed(as.integer(seed.cv) + 7919L)
    }
    n_at <- nrow(X)
    build_folds_at <- function(n_pts, k, strat) {
      if (isTRUE(strat) && n_pts >= 20L) {
        nb <- min(10L, max(2L, floor(n_pts / 5)))
        bins <- cut(as.numeric(y), breaks = nb,
                    include.lowest = TRUE, labels = FALSE)
        fid <- integer(n_pts)
        for (b in seq_len(nb)) {
          idx_b <- which(bins == b)
          if (length(idx_b) == 0L) next
          idx_b <- sample(idx_b)
          fid[idx_b] <- rep_len(seq_len(k), length(idx_b))
        }
        if (any(fid == 0L)) {
          miss <- which(fid == 0L)
          fid[miss] <- rep_len(seq_len(k), length(miss))
        }
        fid
      } else {
        rep_len(seq_len(k), n_pts)[sample.int(n_pts)]
      }
    }

    sigma_grid_sorted <- autotune.grid
    n_sig <- length(sigma_grid_sorted)
    mse_mat_at <- matrix(NA_real_, nrow = n_sig,
                         ncol = nfold_at * ncross_at)
    lam_mat_at <- matrix(NA_real_, nrow = n_sig,
                         ncol = nfold_at * ncross_at)
    lambda_args_nys <- list(
      tol            = 1e-6,
      L0             = .Machine$double.eps,
      L_step         = 10.0,
      U_start_from_n = TRUE
    )
    X_means_loc <- colMeans(X.init)
    X_sds_loc   <- X.init.sd
    nthreads_used_final <- 1L

    col_at <- 1L
    for (cc_at in seq_len(ncross_at)) {
      fid_at <- build_folds_at(n_at, nfold_at, stratify)
      for (k in seq_len(nfold_at)) {
        test_i  <- which(fid_at == k)
        train_i <- which(fid_at != k)
        if (length(test_i) == 0L || length(train_i) < 3L) {
          col_at <- col_at + 1L
          next
        }
        Xs_tr <- Xs[train_i, , drop = FALSE]
        ys_tr <- as.numeric(ys[train_i])
        Xs_te <- Xs[test_i,  , drop = FALSE]
        ys_te <- as.numeric(ys[test_i])
        n_tr  <- nrow(Xs_tr)

        nystrom_m_eff <- if (is.null(nystrom_m)) {
          ceiling(sqrt(n_tr) * 3)
        } else {
          nystrom_m
        }
        fold_landmark_seed <- if (is.null(landmark_seed)) {
          NULL
        } else {
          as.integer(landmark_seed + cc_at * 1000L + k)
        }
        lm_fold <- .resolve_landmarks(landmarks, landmark_method,
                                      nystrom_m_eff, Xs_tr,
                                      X_centers = X_means_loc,
                                      X_scales  = X_sds_loc,
                                      landmark_seed = fold_landmark_seed)
        Z_std <- lm_fold$matrix

        inner <- krls_nystrom_autotune_inner_cpp(
          Xs_tr, Z_std, Xs_te, ys_tr, ys_te,
          sigma_grid_sorted, lambda_args_nys, nystrom_eps,
          nthreads = autotune.nthreads
        )
        mse_mat_at[, col_at] <- as.numeric(inner$mse_per_sigma)
        lam_mat_at[, col_at] <- as.numeric(inner$lambda_per_sigma)
        nthreads_used_final  <- as.integer(inner$nthreads_used)
        col_at <- col_at + 1L
      }
    }

    mean_mse_at <- rowMeans(mse_mat_at, na.rm = TRUE)
    if (all(is.na(mean_mse_at)))
      stop("krls: autotune produced no usable folds.")
    best_si <- which.min(mean_mse_at)
    sd_at  <- apply(mse_mat_at, 1L, function(z) sd(z, na.rm = TRUE))
    nfok   <- rowSums(!is.na(mse_mat_at))
    se_at  <- sd_at / sqrt(pmax(nfok, 1L))
    thresh_at <- mean_mse_at[best_si] + se_at[best_si]
    cand_at   <- which(mean_mse_at <= thresh_at)
    sigma_chosen <- max(sigma_grid_sorted[cand_at])

    ## Refit at chosen sigma on full Xs.
    nystrom_m_full <- if (is.null(nystrom_m)) {
      ceiling(sqrt(nrow(Xs)) * 3)
    } else {
      nystrom_m
    }
    lm_full <- .resolve_landmarks(landmarks, landmark_method,
                                  nystrom_m_full, Xs,
                                  X_centers = X_means_loc,
                                  X_scales  = X_sds_loc,
                                  landmark_seed = landmark_seed)
    nys <- krls_nystrom_fit_cpp(
      Xs, lm_full$matrix, ys, sigma_chosen, lambda_args_nys,
      nystrom_eps, compute_vcov = isTRUE(vcov)
    )
    yfitted_nys <- as.numeric(nys$fitted_std) * as.numeric(y.init.sd) +
                    y.init.mean

    fit <- list(
      coeffs    = nys$coeffs,
      fitted    = yfitted_nys,
      X         = X.init,
      y         = y.init,
      sigma     = sigma_chosen,
      lambda    = nys$lambda,
      X_means   = X_means_loc,
      X_sds     = X_sds_loc,
      y_mean    = y.init.mean,
      y_sd      = as.numeric(y.init.sd),
      approx    = "nystrom",
      landmarks = lm_full$matrix,
      landmark_indices = lm_full$indices,
      landmark_method  = lm_full$method_used,
      nystrom_m        = nys$nystrom_m,
      nystrom_eps      = nystrom_eps,
      W_eigen          = nys$W_eigen,
      Dinvsqrt         = nys$Dinvsqrt,
      Sigma2           = nys$Sigma2,
      vcov_alpha       = if (isTRUE(vcov)) nys$vcov_alpha else NULL,
      autotune         = list(
        grid              = sigma_grid_sorted,
        mse               = mean_mse_at,
        mse_per_sigma     = mean_mse_at,
        winner            = sigma_chosen,
        nfold             = nfold_at,
        ncross            = ncross_at,
        stratify          = isTRUE(stratify),
        seed.cv           = seed.cv,
        mse_per_fold      = mse_mat_at,
        lam_per_fold      = lam_mat_at,
        se_mse            = se_at,
        cv.1se            = TRUE,
        sigma_1se         = sigma_chosen,
        sigma_grid_sorted = sigma_grid_sorted,
        nthreads_used     = nthreads_used_final
      ),
      nystrom_diagnostics = list(
        floored_count = nys$floored_count,
        D_min_raw     = nys$D_min_raw,
        D_max_raw     = nys$D_max_raw
      )
    )
    fit$R2 <- 1 - sum((as.numeric(y.init) - yfitted_nys)^2) /
                  sum((as.numeric(y.init) - mean(as.numeric(y.init)))^2)
    fit$call         <- cl
    fit$na.action    <- na.action
    fit$na.medians   <- na_medians
    fit$factor_info  <- factor_info
    class(fit) <- c("krls_rr", "krls")
    return(fit)
  }

  if (isTRUE(autotune)) {
    nfold_at <- as.integer(nfold)
    if (is.na(nfold_at) || nfold_at <= 0L) nfold_at <- 10L   # Fix 3a
    if (nfold_at > nrow(X))
      stop("krls: autotune nfold (", nfold_at,
           ") exceeds nrow(X) (", nrow(X), ").")
    ## Fix 3b: ncross_at — default 2 when ncross=NULL; respect explicit values.
    ncross_at <- if (is.null(ncross)) 2L else as.integer(ncross)
    if (is.na(ncross_at) || ncross_at < 1L) ncross_at <- 1L

    if (!is.null(seed.cv)) {
      if (exists(".Random.seed", envir = globalenv(),
                 inherits = FALSE)) {
        old_seed_at <- get(".Random.seed", envir = globalenv(),
                            inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed_at,
                       envir = globalenv()), add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()),
                add = TRUE)
      }
      set.seed(as.integer(seed.cv) + 7919L)  # distinct from CV stream
    }
    n_at <- nrow(X)
    ## Helper: build one stratified-or-random fold assignment.
    build_folds_at <- function(n_pts, k, strat) {
      if (isTRUE(strat) && n_pts >= 20L) {
        nb <- min(10L, max(2L, floor(n_pts / 5)))
        bins <- cut(as.numeric(y), breaks = nb,
                    include.lowest = TRUE, labels = FALSE)
        fid <- integer(n_pts)
        for (b in seq_len(nb)) {
          idx_b <- which(bins == b)
          if (length(idx_b) == 0L) next
          idx_b <- sample(idx_b)
          fid[idx_b] <- rep_len(seq_len(k), length(idx_b))
        }
        if (any(fid == 0L)) {
          miss <- which(fid == 0L)
          fid[miss] <- rep_len(seq_len(k), length(miss))
        }
        fid
      } else {
        rep_len(seq_len(k), n_pts)[sample.int(n_pts)]
      }
    }
    ## MSE matrix: rows = sigma candidates, cols = nfold_at * ncross_at.
    n_sig <- length(autotune.grid)
    sigma_grid_sorted <- autotune.grid   # already sorted+unique above
    mse_mat_at <- matrix(NA_real_, nrow = n_sig,
                         ncol = nfold_at * ncross_at)
    lam_mat_at <- matrix(NA_real_, nrow = n_sig,
                         ncol = nfold_at * ncross_at)
    ## Use C++ inner only on the fast unweighted path with default
    ## L/U/tol (the inner uses L0=eps + L_step=10 + U=n which matches
    ## .krls_lambdasearch defaults exactly; non-default L/U/tol or
    ## non-null weights fall back to the per-cell R loop for parity).
    use_inner_cpp <- is.null(sqw) && is.null(L) && is.null(U) &&
                       is.null(tol)
    lambda_args <- list(
      tol            = 1e-6,
      L0             = .Machine$double.eps,
      L_step         = 10.0,
      U_start_from_n = TRUE
    )
    nthreads_used_final <- 1L
    col_at <- 1L
    for (cc_at in seq_len(ncross_at)) {
      ## Fresh fold partition for each repetition.
      fid_at <- build_folds_at(n_at, nfold_at, stratify)
      for (k in seq_len(nfold_at)) {
        test_i  <- which(fid_at == k)
        train_i <- which(fid_at != k)
        if (length(test_i) == 0L || length(train_i) < 3L) {
          col_at <- col_at + 1L
          next
        }
        Xs_tr <- Xs[train_i, , drop = FALSE]
        ys_tr <- ys[train_i]
        Xs_te <- Xs[test_i,  , drop = FALSE]
        ys_te <- ys[test_i]
        if (use_inner_cpp) {
          ## Shared squared-distance matrices for this fold; the C++
          ## inner sweeps the sigma grid in parallel reusing them.
          D_tr <- krls_pairwise_sqdist_cpp(Xs_tr, Xs_tr)
          D_te <- krls_pairwise_sqdist_cpp(Xs_te, Xs_tr)
          inner <- krls_autotune_inner_cpp(
            D_tr, D_te, ys_tr, ys_te,
            sigma_grid_sorted, lambda_args,
            nthreads = autotune.nthreads
          )
          mse_mat_at[, col_at] <- as.numeric(inner$mse_per_sigma)
          lam_mat_at[, col_at] <- as.numeric(inner$lambda_per_sigma)
          nthreads_used_final  <- as.integer(inner$nthreads_used)
        } else {
          for (si in seq_along(autotune.grid)) {
            sig_s <- autotune.grid[si]
            K_tr <- krls_kernel_cpp(Xs_tr, sig_s)
            K_te <- krls_kernel_pred_cpp(Xs_te, Xs_tr, sig_s)
            if (!is.null(sqw)) {
              sqw_tr <- sqrt(w_norm[train_i])
              K_tr_use  <- K_tr * tcrossprod(sqw_tr)
              ys_tr_use <- sqw_tr * ys_tr
            } else {
              sqw_tr <- NULL
              K_tr_use  <- K_tr
              ys_tr_use <- ys_tr
            }
            eo_k    <- krls_eig_cpp(K_tr_use)
            dvals_k <- as.numeric(eo_k$values)
            V_k     <- eo_k$vectors
            Vsq_k   <- krls_vsq_cpp(V_k)
            Vty_k   <- as.numeric(crossprod(V_k, ys_tr_use))
            ## Per-fold LOO closed-form lambda search for this sigma cell.
            lam_k <- tryCatch(
              .krls_lambdasearch(dvals_k, V_k, Vsq_k, Vty_k,
                                 n_y = length(train_i),
                                 L = L, U = U, tol = tol, noisy = FALSE),
              error = function(e) NA_real_)
            if (is.na(lam_k)) next
            sol_k   <- krls_solve_cpp(dvals_k, V_k, Vsq_k, Vty_k, lam_k)
            c_solve <- as.numeric(sol_k$coeffs)
            c_fold  <- if (!is.null(sqw_tr)) sqw_tr * c_solve else c_solve
            yhat_te <- as.numeric(K_te %*% c_fold)
            if (!is.null(sqw)) {
              w_te <- w_norm[test_i]
              mse_mat_at[si, col_at] <- sum(w_te * (ys_te - yhat_te)^2) /
                max(sum(w_te), .Machine$double.eps)
            } else {
              mse_mat_at[si, col_at] <- mean((ys_te - yhat_te)^2)
            }
            lam_mat_at[si, col_at] <- lam_k
          }
        }
        col_at <- col_at + 1L
      }
    }
    mean_mse_at <- rowMeans(mse_mat_at, na.rm = TRUE)
    if (all(is.na(mean_mse_at)))
      stop("krls: autotune produced no usable folds.")
    best_si <- which.min(mean_mse_at)
    ## Fix 3c: 1-SE rule — select the LARGEST sigma within 1 SE of min.
    ## Wider kernel = more smoothing = less overfit.
    sd_at  <- apply(mse_mat_at, 1L, function(z) sd(z, na.rm = TRUE))
    nfok   <- rowSums(!is.na(mse_mat_at))
    se_at  <- sd_at / sqrt(pmax(nfok, 1L))
    thresh_at <- mean_mse_at[best_si] + se_at[best_si]
    cand_at   <- which(mean_mse_at <= thresh_at)
    sigma <- max(autotune.grid[cand_at])   # largest sigma in the 1-SE band
    autotune_info <- list(
      grid              = autotune.grid,
      mse               = mean_mse_at,
      mse_per_sigma     = mean_mse_at,
      winner            = sigma,
      nfold             = nfold_at,
      ncross            = ncross_at,
      stratify          = isTRUE(stratify),
      seed.cv           = seed.cv,
      mse_per_fold      = mse_mat_at,
      lam_per_fold      = lam_mat_at,
      se_mse            = se_at,
      cv.1se            = TRUE,
      sigma_1se         = sigma,
      sigma_grid_sorted = sigma_grid_sorted,
      nthreads_used     = nthreads_used_final
    )
    if (trace > 1) {
      cat("Autotune sigma winner (1-SE):", round(sigma, 4),
          " (mse=", round(mean_mse_at[best_si], 6), ")\n", sep = "")
    }
  }

  ## --- Nystrom non-autotune routing (Phase 2 T4) --------------------
  if (identical(approx, "nystrom") && !isTRUE(autotune)) {
    if (!is.null(w_norm)) {
      stop("krls: approx = 'nystrom' does not yet support `weights`.",
           call. = FALSE)
    }
    X_means_loc <- colMeans(X.init)
    X_sds_loc   <- X.init.sd
    lm <- .resolve_landmarks(landmarks, landmark_method, nystrom_m,
                             Xs, X_centers = X_means_loc,
                             X_scales = X_sds_loc,
                             landmark_seed = landmark_seed)
    Z_std <- lm$matrix

    lambda_args_nys <- list(
      tol = 1e-6, L0 = .Machine$double.eps,
      L_step = 10.0, U_start_from_n = TRUE
    )

    nys <- krls_nystrom_fit_cpp(
      Xs, Z_std, ys, sigma, lambda_args_nys,
      nystrom_eps, compute_vcov = isTRUE(vcov)
    )

    yfitted_nys <- as.numeric(nys$fitted_std) * as.numeric(y.init.sd) +
                     y.init.mean

    fit <- list(
      coeffs    = nys$coeffs,
      fitted    = yfitted_nys,
      X         = X.init,
      y         = y.init,
      sigma     = sigma,
      lambda    = nys$lambda,
      X_means   = X_means_loc,
      X_sds     = X_sds_loc,
      y_mean    = y.init.mean,
      y_sd      = as.numeric(y.init.sd),
      approx    = "nystrom",
      landmarks = Z_std,
      landmark_indices = lm$indices,
      landmark_method  = lm$method_used,
      nystrom_m        = nys$nystrom_m,
      nystrom_eps      = nystrom_eps,
      W_eigen          = nys$W_eigen,
      Dinvsqrt         = nys$Dinvsqrt,
      Sigma2           = nys$Sigma2,
      vcov_alpha       = if (isTRUE(vcov)) nys$vcov_alpha else NULL,
      nystrom_diagnostics = list(
        floored_count = nys$floored_count,
        D_min_raw     = nys$D_min_raw,
        D_max_raw     = nys$D_max_raw
      )
    )

    fit$R2 <- 1 - sum((as.numeric(y.init) - yfitted_nys)^2) /
                  sum((as.numeric(y.init) - mean(as.numeric(y.init)))^2)

    fit$call         <- cl
    fit$na.action    <- na.action
    fit$na.medians   <- na_medians
    fit$factor_info  <- factor_info

    class(fit) <- c("krls_rr", "krls")
    return(fit)
  }

  ## --- kernel + eigendecomposition ---------------------------------
  K <- krls_kernel_cpp(Xs, sigma)
  if (!is.null(sqw)) {
    # Twist the kernel: K_w = D K D
    K_w <- K * tcrossprod(sqw)
    eo <- krls_eig_cpp(K_w)
    dvals <- as.numeric(eo$values)
    V     <- eo$vectors
    Vsq   <- krls_vsq_cpp(V)
    ys_w  <- sqw * ys
    Vty   <- as.numeric(crossprod(V, ys_w))
    yty   <- sum(ys_w^2)
  } else {
    eo <- krls_eig_cpp(K)
    dvals <- as.numeric(eo$values)
    V     <- eo$vectors
    Vsq   <- krls_vsq_cpp(V)
    Vty   <- as.numeric(crossprod(V, ys))
    yty   <- sum(ys^2)
  }

  ## --- eigentruncation handling (optional) -------------------------
  if (!is.null(eigtrunc)) {
    lastkeeper <- max(which(dvals >= eigtrunc * dvals[1L]))
    lastkeeper <- max(1L, lastkeeper)
    if (lastkeeper < length(dvals)) {
      dvals <- dvals[seq_len(lastkeeper)]
      V     <- V[, seq_len(lastkeeper), drop = FALSE]
      Vsq   <- Vsq[, seq_len(lastkeeper), drop = FALSE]
      Vty   <- Vty[seq_len(lastkeeper)]
    }
  }

  ## --- lambda selection --------------------------------------------
  # Phase 3: lambda.method = c('loo', 'gcv', 'cv'). LOO is the
  # back-compat default and uses the closed-form Hainmueller-Hazlett
  # identity. GCV uses the closed-form Craven-Wahba criterion on the
  # same eigenbasis. CV does K-fold (or repeated K-fold) over an
  # explicit lambda.grid and picks argmin held-out MSE (or 1-SE),
  # optionally seeded for reproducibility.
  cv_diag <- NULL
  if (is.null(lambda)) {
    if (lambda.method %in% c("loo", "gcv")) {
      lambda <- .krls_lambdasearch(dvals = dvals, V = V, Vsq = Vsq,
                                   Vty = Vty, n_y = nrow(y),
                                   yty = yty,
                                   L = L, U = U, tol = tol,
                                   method = lambda.method,
                                   noisy = isTRUE(trace > 2))
      if (trace > 1) {
        label <- if (lambda.method == "loo") "Loo-Loss" else "GCV"
        cat("Lambda that minimizes ", label, " is:",
            round(lambda, 5), "\n", sep = "")
      }
    } else {
      # CV path. Default nfold = 5 if user asked lambda.method='cv' but
      # didn't specify nfold. nfold = nrow(X) reproduces LOO via CV.
      nfold_int <- as.integer(nfold)
      if (is.na(nfold_int) || nfold_int <= 0L) nfold_int <- 5L
      if (nfold_int > nrow(X))
        stop("krls: nfold (", nfold_int, ") exceeds nrow(X) (",
             nrow(X), ").")
      ncross_int <- if (is.null(ncross)) 1L else as.integer(ncross)
      if (is.na(ncross_int) || ncross_int < 1L) ncross_int <- 1L

      # Default lambda.grid: 30 log-spaced points within the standard
      # (L, U) bracket from the closed-form search (which is itself
      # derived from the eigenvalues so the grid is data-adaptive).
      if (is.null(lambda.grid)) {
        # Reuse the same L/U bracket the golden search would use.
        U_grid <- if (is.null(U)) {
          Uloc <- nrow(y)
          while (sum(dvals / (dvals + Uloc)) < 1) Uloc <- Uloc - 1
          Uloc
        } else U
        L_grid <- if (is.null(L)) {
          q  <- which.min(abs(dvals - (max(dvals) / 1000)))
          Lloc  <- .Machine$double.eps
          while (sum(dvals / (dvals + Lloc)) > q) Lloc <- Lloc * 10
          Lloc
        } else max(L, .Machine$double.eps)
        if (U_grid <= L_grid) {
          U_grid <- max(L_grid * 10, L_grid + 1)
        }
        lambda.grid <- exp(seq(log(max(L_grid, .Machine$double.eps)),
                               log(U_grid), length.out = 30L))
      } else {
        stopifnot(is.numeric(lambda.grid),
                  length(lambda.grid) >= 1L,
                  all(lambda.grid > 0))
        lambda.grid <- sort(unique(as.numeric(lambda.grid)))
      }

      # RNG protection so CV only consumes the user's RNG when
      # seed.cv = NULL.
      if (!is.null(seed.cv)) {
        if (exists(".Random.seed", envir = globalenv(),
                   inherits = FALSE)) {
          old_seed <- get(".Random.seed", envir = globalenv(),
                          inherits = FALSE)
          on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                  add = TRUE)
        } else {
          on.exit(rm(list = ".Random.seed", envir = globalenv()),
                  add = TRUE)
        }
        set.seed(as.integer(seed.cv))
      }

      n_cv <- nrow(X)
      # Stratification for continuous y: use quantile bins (10 bins or
      # fewer if n < 50).  For a small-y, fall back to plain shuffle.
      build_folds <- function(n_pts, k, strat) {
        if (isTRUE(strat) && n_pts >= 20L) {
          nb <- min(10L, max(2L, floor(n_pts / 5)))
          bins <- cut(as.numeric(y), breaks = nb, include.lowest = TRUE,
                      labels = FALSE)
          # Distribute each bin's indices round-robin across folds.
          fid <- integer(n_pts)
          for (b in seq_len(nb)) {
            idx_b <- which(bins == b)
            if (length(idx_b) == 0L) next
            idx_b <- sample(idx_b)  # shuffle within bin
            fid[idx_b] <- rep_len(seq_len(k), length(idx_b))
          }
          # Any unassigned (NA bin) get round-robin too.
          if (any(fid == 0L)) {
            miss <- which(fid == 0L)
            fid[miss] <- rep_len(seq_len(k), length(miss))
          }
          fid
        } else {
          rep_len(seq_len(k), n_pts)[sample.int(n_pts)]
        }
      }

      mse_mat <- matrix(NA_real_, nrow = length(lambda.grid),
                        ncol = nfold_int * ncross_int)
      col_idx <- 1L
      for (cc in seq_len(ncross_int)) {
        fid <- build_folds(n_cv, nfold_int, stratify)
        for (k in seq_len(nfold_int)) {
          test_i  <- which(fid == k)
          train_i <- which(fid != k)
          if (length(test_i) == 0L || length(train_i) < 3L) next
          Xs_tr <- Xs[train_i, , drop = FALSE]
          ys_tr <- ys[train_i]
          Xs_te <- Xs[test_i,  , drop = FALSE]
          ys_te <- ys[test_i]
          # Fold-local kernel + eigendecomp (the shared piece across
          # all lambda candidates per fold).
          K_tr <- krls_kernel_cpp(Xs_tr, sigma)
          K_te <- krls_kernel_pred_cpp(Xs_te, Xs_tr, sigma)
          if (!is.null(w_norm)) {
            sqw_tr <- sqrt(w_norm[train_i])
            K_tr_use <- K_tr * tcrossprod(sqw_tr)
            ys_tr_use <- sqw_tr * ys_tr
          } else {
            sqw_tr <- NULL
            K_tr_use <- K_tr
            ys_tr_use <- ys_tr
          }
          eo_k <- krls_eig_cpp(K_tr_use)
          dvals_k <- as.numeric(eo_k$values)
          V_k     <- eo_k$vectors
          Vsq_k   <- krls_vsq_cpp(V_k)
          Vty_k   <- as.numeric(crossprod(V_k, ys_tr_use))
          for (li in seq_along(lambda.grid)) {
            lam_l <- lambda.grid[li]
            sol_k <- krls_solve_cpp(dvals_k, V_k, Vsq_k, Vty_k, lam_l)
            c_solve <- as.numeric(sol_k$coeffs)
            c_fold  <- if (!is.null(sqw_tr)) sqw_tr * c_solve else c_solve
            yhat_te <- as.numeric(K_te %*% c_fold)
            err_v   <- (ys_te - yhat_te)
            if (!is.null(w_norm)) {
              w_te <- w_norm[test_i]
              mse_mat[li, col_idx] <-
                sum(w_te * err_v^2) / max(sum(w_te), .Machine$double.eps)
            } else {
              mse_mat[li, col_idx] <- mean(err_v^2)
            }
          }
          col_idx <- col_idx + 1L
        }
      }
      # Average across folds / repeats. Use min().
      mean_mse <- rowMeans(mse_mat, na.rm = TRUE)
      if (all(is.na(mean_mse)))
        stop("krls: lambda.method='cv' produced no usable folds; ",
             "check nfold / stratify.")
      best_i <- which.min(mean_mse)
      lambda <- lambda.grid[best_i]
      if (isTRUE(cv.1se)) {
        # 1-SE rule: pick the largest lambda whose mean is within 1 SE
        # of the best. SE = sd across folds / sqrt(n_folds).
        sd_mse <- apply(mse_mat, 1L, function(z) sd(z, na.rm = TRUE))
        nfok <- rowSums(!is.na(mse_mat))
        se_mse <- sd_mse / sqrt(pmax(nfok, 1L))
        thresh <- mean_mse[best_i] + se_mse[best_i]
        cand <- which(mean_mse <= thresh)
        if (length(cand)) lambda <- max(lambda.grid[cand])
      }
      cv_diag <- list(
        method      = "cv",
        nfold       = nfold_int,
        ncross      = ncross_int,
        stratify    = isTRUE(stratify),
        seed.cv     = seed.cv,
        lambda.grid = lambda.grid,
        mean_mse    = mean_mse,
        mse_per_fold = mse_mat,
        best_idx    = best_i,
        cv.1se      = isTRUE(cv.1se)
      )
      if (trace > 1) {
        cat("CV lambda (", nfold_int, "-fold x ", ncross_int, " cross):",
            round(lambda, 6), "\n", sep = "")
      }
    }
  } else {
    stopifnot(is.numeric(lambda), length(lambda) == 1, lambda > 0)
  }

  ## --- solve for c, LOOe, diag(G^{-1}) -----------------------------
  sol     <- krls_solve_cpp(dvals, V, Vsq, Vty, lambda)
  coeffs_solve <- as.numeric(sol$coeffs)
  Le      <- sol$Le
  if (!is.null(sqw)) {
    # Un-twist coefficients: c_final = D c_solve so K c_final has the
    # right meaning on the original (un-D-scaled) input.
    coeffs <- sqw * coeffs_solve
    yfit_s <- as.numeric(K %*% coeffs)
  } else {
    coeffs <- coeffs_solve
    yfit_s  <- as.numeric(K %*% coeffs)
  }

  ## --- coefficient covariance --------------------------------------
  vcov.c <- NULL
  vcov.fitted <- NULL
  if (!is.null(sqw)) {
    # Weighted residual variance: sum w * (ys - yfit_s)^2 / sum(w).
    sigma2 <- as.numeric(sum(w_norm * (ys - yfit_s)^2) / sum(w_norm))
  } else {
    sigma2 <- as.numeric((1 / n) * crossprod(ys - yfit_s))
  }
  if (isTRUE(vcov)) {
    w <- sigma2 / (dvals + lambda)^2
    vcovmatc_solve <- V %*% (t(V) * w)
    if (!is.null(sqw)) {
      # Un-twist: vcov(c_final) = D %*% vcov(c_solve) %*% D
      vcovmatc <- vcovmatc_solve * tcrossprod(sqw)
    } else {
      vcovmatc <- vcovmatc_solve
    }
    vcovmatyhat <- crossprod(K, vcovmatc %*% K)
  }

  ## --- derivatives --------------------------------------------------
  derivmat <- avgderiv <- varavgderivmat <- NULL
  if (isTRUE(derivative)) {
    derivmat_s <- krls_deriv_cpp(Xs, K, coeffs, sigma)
    if (!is.null(sqw)) {
      avgderiv_s <- matrix(
        colSums(derivmat_s * w_norm) / sum(w_norm),
        nrow = 1L)
    } else {
      avgderiv_s <- matrix(colMeans(derivmat_s), nrow = 1L)
    }
    varavg_s   <- krls_avg_deriv_var_cpp(Xs, K, V, dvals, sigma, lambda, sigma2)
    colnames(derivmat_s)  <- colnames(X)
    colnames(avgderiv_s)  <- colnames(X)
    scale_vec <- as.numeric(y.init.sd) / X.init.sd
    derivmat  <- sweep(derivmat_s, 2L, scale_vec, "*")
    if (!is.null(sqw)) {
      avgderiv  <- matrix(
        colSums(derivmat * w_norm) / sum(w_norm),
        nrow = 1L,
        dimnames = list(NULL, colnames(X)))
    } else {
      avgderiv  <- matrix(colMeans(derivmat), nrow = 1L,
                          dimnames = list(NULL, colnames(X)))
    }
    varavgderivmat <- matrix(scale_vec^2 * as.numeric(varavg_s),
                             nrow = 1L,
                             dimnames = list(NULL, colnames(X)))
  }

  ## --- unstandardise fitted, vcov ----------------------------------
  yfitted <- yfit_s * as.numeric(y.init.sd) + y.init.mean
  if (isTRUE(vcov)) {
    vcov.c      <- (as.numeric(y.init.sd)^2) * vcovmatc
    vcov.fitted <- (as.numeric(y.init.sd)^2) * vcovmatyhat
  }

  Looe <- Le * as.numeric(y.init.sd)
  R2   <- 1 - (var(as.numeric(y.init - yfitted)) / (as.numeric(y.init.sd)^2))

  # --- residual variance model (Phase 6) ---------------------------
  # Effective degrees of freedom for the smoother: trace of the hat
  # matrix H = K (K + lambda I)^{-1} = V diag(d/(d+lambda)) V'.
  # tr(H) = sum_i d_i / (d_i + lambda).
  varmod_info <- NULL
  if (varmod != "none") {
    effdf <- sum(dvals / (dvals + lambda))
    resid_raw <- as.numeric(y.init) - yfitted
    if (!is.null(w_norm)) {
      rss_w <- sum(w_norm * resid_raw^2)
      df_eff <- max(sum(w_norm) - effdf, 1)
      sigma_hat <- sqrt(rss_w / df_eff)
    } else {
      rss <- sum(resid_raw^2)
      df_eff <- max(n - effdf, 1)
      sigma_hat <- sqrt(rss / df_eff)
    }
    varmod_info <- list(
      type      = varmod,
      sigma_hat = sigma_hat,
      effdf     = effdf,
      df        = df_eff
    )
  }

  binaryindicator <- matrix(FALSE, 1L, d,
                            dimnames = list(NULL, colnames(X)))

  out <- list(K = K, coeffs = coeffs, Looe = Looe, fitted = yfitted,
              X = X.init, y = y.init,
              sigma = sigma, lambda = lambda, R2 = R2,
              derivatives = derivmat,
              avgderivatives = avgderiv,
              var.avgderivatives = varavgderivmat,
              vcov.c = vcov.c, vcov.fitted = vcov.fitted,
              binaryindicator = binaryindicator)

  # Bookkeeping for formula method / predict() / update()
  out$call <- cl
  out$na.action <- na.action
  out$na.medians <- na_medians
  out$factor_info <- factor_info
  out$weights <- if (!is.null(w_norm)) w_norm else NULL
  if (!is.null(cv_diag)) out$cv <- cv_diag
  if (!is.null(autotune_info)) out$autotune <- autotune_info
  if (!is.null(varmod_info)) out$varmod <- varmod_info

  class(out) <- c("krls_rr", "krls")

  if (isTRUE(derivative) && isTRUE(binary)) {
    out <- .krls_fd_binary(out)
  }

  # --- bagging (Phase 5) -------------------------------------------
  # Bootstrap n_boot row samples; for each, refit krls inheriting the
  # central fit's sigma + lambda (so we don't re-run autotune / CV /
  # LOO search per replicate). Store the minimum predict-ready state
  # per replicate.
  n_boot <- as.integer(n.boot)
  if (is.na(n_boot) || n_boot < 0L) n_boot <- 0L
  if (n_boot > 0L) {
    sig_b   <- out$sigma
    lam_b   <- out$lambda
    n_full  <- nrow(X)
    if (!is.null(seed.cv)) {
      if (exists(".Random.seed", envir = globalenv(),
                 inherits = FALSE)) {
        old_seed_b <- get(".Random.seed", envir = globalenv(),
                           inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed_b,
                       envir = globalenv()), add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()),
                add = TRUE)
      }
      set.seed(as.integer(seed.cv) + 31337L)
    }
    boot_reps <- vector("list", n_boot)
    idx_list  <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n_full, n_full, replace = TRUE)
      idx_list[[b]] <- idx
      X_b <- X[idx, , drop = FALSE]
      y_b <- as.numeric(y)[idx]
      w_b <- if (!is.null(w_norm)) {
        wb <- w_norm[idx]
        wb / mean(wb)
      } else NULL
      # Fit on the bootstrap sample with the inherited (sigma, lambda).
      # Derivatives / vcov are dropped to keep memory bounded.
      f_b <- krls.default(
        X = X_b, y = y_b,
        sigma = sig_b, lambda = lam_b,
        derivative = FALSE, binary = FALSE, vcov = FALSE,
        weights = w_b,
        lambda.method = "loo",
        autotune = FALSE,
        varmod = "none",
        n.boot = 0L,
        na.action = na.action,
        trace = 0L, nthreads = nthreads
      )
      # Strip heavy fields; keep what predict needs:
      #   - X (training rows in raw scale, used to rebuild Xs_b in predict),
      #   - y (for mean/sd),
      #   - coeffs, sigma, lambda.
      boot_reps[[b]] <- list(
        X      = f_b$X,
        y      = f_b$y,
        coeffs = f_b$coeffs,
        sigma  = f_b$sigma,
        lambda = f_b$lambda
      )
    }
    out$boot <- list(replicates = boot_reps, n.boot = n_boot,
                     idx = idx_list)
  }

  if (trace > 0 && isTRUE(derivative)) {
    av <- setNames(as.vector(out$avgderivatives),
                   colnames(out$avgderivatives))
    cat("\n Average Marginal Effects:\n \n"); print(av)
    cat("\n Quartiles of Marginal Effects:\n \n")
    print(apply(out$derivatives, 2L, quantile,
                probs = c(0.25, 0.5, 0.75)))
  }
  out
}

#' @rdname krls
#' @param object A fitted `"krls"` object.
#' @param newdata A numeric matrix or data frame with the same columns
#'   as the training `X` (or matching the training formula).  `NULL`
#'   returns `object$fitted`.
#' @param se.fit Logical.  If `TRUE`, return pointwise standard errors
#'   of the predictions.  Requires the fit was created with
#'   `vcov = TRUE`.
#' @param interval Prediction-interval mode. `"none"` (default)
#'   returns point predictions only; `"pint"` returns lower/upper
#'   bounds at confidence level `level`, requires the fit was
#'   created with `varmod = "const"` (or other non-`"none"`
#'   `varmod`).
#' @param level Confidence level for `interval = "pint"`. Default
#'   `0.95`.
#' @export
predict.krls_rr <- function(object, newdata = NULL, se.fit = FALSE,
                              interval = c("none", "pint"),
                              level = 0.95, ...) {
  interval <- match.arg(interval)
  if (!inherits(object, "krls")) {
    stop("object is not of class 'krls'")
  }
  if (isTRUE(se.fit) && is.null(object$vcov.c)) {
    stop("refit with krls(..., vcov = TRUE) to compute standard errors")
  }
  if (interval == "pint") {
    if (is.null(object$varmod) || is.null(object$varmod$sigma_hat))
      stop("krls: interval='pint' requires the fit was built with",
           " varmod = 'const'. Refit krls() with varmod = 'const' and",
           " try again.", call. = FALSE)
    if (is.null(object$vcov.c))
      stop("krls: interval='pint' requires vcov = TRUE at fit time.",
           call. = FALSE)
  }

  # --- Nystrom predict branch (Phase 2 T4) ---------------------------
  if (identical(object$approx, "nystrom")) {
    if (is.null(newdata)) {
      return(list(fit = matrix(object$fitted, ncol = 1L),
                  se.fit = NULL, vcov.fit = NULL,
                  newdata = NULL, newdataK = NULL))
    }
    xnew <- .krls_build_design(object, newdata)
    storage.mode(xnew) <- "double"
    if (ncol(object$X) != ncol(xnew)) {
      stop("ncol(newdata) (", ncol(xnew),
           ") differs from ncol(X) from fitted krls object (",
           ncol(object$X), ")")
    }
    if (any(is.na(xnew))) {
      if (!is.null(object$na.medians)) {
        for (j in seq_len(ncol(xnew))) {
          na_j <- is.na(xnew[, j])
          if (any(na_j)) xnew[na_j, j] <- object$na.medians[j]
        }
        warning("krls: median-imputed ", sum(is.na(xnew)),
                " missing value(s) in newdata using stored training medians.",
                call. = FALSE)
      } else {
        stop("krls: newdata has NA values but the fit has no stored medians",
             " (training used na.action='omit'). Drop or impute the rows",
             " yourself before calling predict().")
      }
    }
    Xn_std <- scale(xnew, center = object$X_means, scale = object$X_sds)
    attr(Xn_std, "scaled:center") <- NULL
    attr(Xn_std, "scaled:scale")  <- NULL
    Xn_std <- matrix(as.numeric(Xn_std), nrow(xnew), ncol(xnew))
    yhat_std <- as.numeric(krls_nystrom_predict_cpp(
      Xn_std, object$landmarks, as.numeric(object$coeffs), object$sigma
    ))
    yhat <- yhat_std * object$y_sd + object$y_mean
    return(list(fit = matrix(yhat, ncol = 1L),
                se.fit = NULL, vcov.fit = NULL,
                newdata = Xn_std, newdataK = NULL))
  }

  if (is.null(newdata)) {
    # For the NULL newdata fast-path, build the PI from $fitted +
    # the stored vcov.fitted (on training points).
    if (interval == "pint") {
      vh <- object$vcov.fitted
      se_t <- if (is.null(vh)) NA_real_ else sqrt(diag(vh))
      sh <- object$varmod$sigma_hat
      df <- object$varmod$df
      tq <- qt(1 - (1 - level) / 2, df = max(df, 1))
      half <- tq * sqrt(se_t^2 + sh^2)
      pi_mat <- cbind(fit = as.numeric(object$fitted),
                      lwr = as.numeric(object$fitted) - half,
                      upr = as.numeric(object$fitted) + half)
      return(pi_mat)
    }
    return(list(fit = matrix(object$fitted, ncol = 1L),
                se.fit = NULL, vcov.fit = NULL,
                newdata = NULL, newdataK = NULL))
  }

  # Build the design matrix to score against.
  xnew <- .krls_build_design(object, newdata)

  storage.mode(xnew) <- "double"
  if (ncol(object$X) != ncol(xnew)) {
    stop("ncol(newdata) (", ncol(xnew),
         ") differs from ncol(X) from fitted krls object (",
         ncol(object$X), ")")
  }

  # Re-impute NAs in newdata using stored medians (matrix path only;
  # the data-frame / formula path imputed inside .krls_build_design).
  if (any(is.na(xnew))) {
    if (!is.null(object$na.medians)) {
      for (j in seq_len(ncol(xnew))) {
        na_j <- is.na(xnew[, j])
        if (any(na_j)) xnew[na_j, j] <- object$na.medians[j]
      }
      warning("krls: median-imputed ", sum(is.na(xnew)),
              " missing value(s) in newdata using stored training medians.",
              call. = FALSE)
    } else {
      stop("krls: newdata has NA values but the fit has no stored medians",
           " (training used na.action='omit'). Drop or impute the rows",
           " yourself before calling predict().")
    }
  }

  # Column-name handling for the bare-matrix path (only applies when
  # no factor_info / terms drove the rebuild above).
  if (is.null(object$factor_info) && is.null(object$terms)) {
    cn_train <- colnames(object$X)
    cn_new   <- colnames(xnew)
    if (!is.null(cn_train) && !is.null(cn_new)) {
      if (!setequal(cn_train, cn_new)) {
        stop("colnames(newdata) do not match colnames(X) from fit: ",
             "training = c(", paste(shQuote(cn_train), collapse = ", "),
             "), newdata = c(", paste(shQuote(cn_new), collapse = ", "),
             ")")
      }
      if (!identical(cn_train, cn_new)) {
        xnew <- xnew[, cn_train, drop = FALSE]
      }
    }
  }

  if (nrow(object$X) < 2L) {
    stop("object$X has fewer than 2 rows; sd cannot be recomputed")
  }
  Xmeans <- colMeans(object$X)
  Xsd    <- apply(object$X, 2L, sd)
  if (any(!is.finite(Xsd)) || any(Xsd == 0)) {
    stop("at least one stored X column has zero or non-finite sd; ",
         "object$X may have been mutated after fit")
  }
  Xs     <- scale(object$X, center = Xmeans, scale = Xsd)
  Xs     <- matrix(Xs, nrow(object$X), ncol(object$X))
  Xn     <- scale(xnew, center = Xmeans, scale = Xsd)
  Xn     <- matrix(Xn, nrow(xnew), ncol(xnew))

  Knew   <- krls_kernel_pred_cpp(Xn, Xs, object$sigma)
  yhat   <- Knew %*% object$coeffs

  vcov.fit <- se.fit.out <- NULL
  need_se <- isTRUE(se.fit) || interval == "pint"
  if (need_se) {
    vcov.c.raw  <- object$vcov.c * as.vector(1 / var(as.vector(object$y)))
    vcov.fitted <- tcrossprod(Knew %*% vcov.c.raw, Knew)
    sd_y2       <- as.numeric(apply(object$y, 2L, sd))^2
    vcov.fit    <- sd_y2 * vcov.fitted
    se.fit.out  <- matrix(sqrt(diag(vcov.fit)), ncol = 1L)
  }
  sd_y_train <- as.numeric(apply(object$y, 2L, sd))
  mean_y_train <- mean(object$y)
  yhat <- (yhat * sd_y_train) + mean_y_train

  # --- Bagging mean (Phase 5) --------------------------------------
  # If the fit carries bootstrap replicates, average the central
  # prediction with the per-replicate predictions on the same newdata.
  bag_sd <- NULL
  if (!is.null(object$boot) && !is.null(object$boot$replicates)) {
    reps <- object$boot$replicates
    nrep <- length(reps)
    if (nrep > 0L) {
      preds <- matrix(NA_real_, nrow = length(yhat), ncol = nrep + 1L)
      preds[, 1L] <- as.numeric(yhat)
      for (b in seq_along(reps)) {
        Rb <- reps[[b]]
        Xmeans_b <- colMeans(Rb$X)
        Xsd_b    <- apply(Rb$X, 2L, sd)
        # Skip replicates with degenerate X (zero sd).
        if (any(!is.finite(Xsd_b)) || any(Xsd_b == 0)) next
        Xs_b <- scale(Rb$X, center = Xmeans_b, scale = Xsd_b)
        Xs_b <- matrix(Xs_b, nrow(Rb$X), ncol(Rb$X))
        Xn_b <- scale(xnew, center = Xmeans_b, scale = Xsd_b)
        Xn_b <- matrix(Xn_b, nrow(xnew), ncol(xnew))
        Kn_b <- krls_kernel_pred_cpp(Xn_b, Xs_b, Rb$sigma)
        yhat_b <- as.numeric(Kn_b %*% Rb$coeffs)
        sd_yb  <- as.numeric(apply(Rb$y, 2L, sd))
        yhat_b <- yhat_b * sd_yb + mean(Rb$y)
        preds[, b + 1L] <- yhat_b
      }
      bag_mean <- rowMeans(preds, na.rm = TRUE)
      bag_sd   <- apply(preds, 1L,
                        function(z) sd(z, na.rm = TRUE))
      yhat <- matrix(bag_mean, ncol = 1L)
    }
  }

  if (interval == "pint") {
    sh <- object$varmod$sigma_hat
    df <- object$varmod$df
    tq <- qt(1 - (1 - level) / 2, df = max(df, 1))
    se_t <- as.numeric(se.fit.out)
    half <- tq * sqrt(se_t^2 + sh^2)
    pi_mat <- cbind(fit = as.numeric(yhat),
                    lwr = as.numeric(yhat) - half,
                    upr = as.numeric(yhat) + half)
    return(pi_mat)
  }
  list(fit = yhat, se.fit = se.fit.out, vcov.fit = vcov.fit,
       newdata = Xn, newdataK = Knew,
       bag_sd = bag_sd)
}

## -----------------------------------------------------------------
## Helpers (not exported).
## -----------------------------------------------------------------

# Build a numeric design matrix from `newdata`, honouring the fit's
# stored factor / terms / column metadata. Ports BUG-004 (OOV factor),
# BUG-010 (factor NA), and BUG-012 (derived terms) from predict.ares.
.krls_build_design <- function(object, newdata) {
  if (is.data.frame(newdata)) {
    if (!is.null(object$factor_info)) {
      fi <- object$factor_info
      if (!all(fi$orig_names %in% colnames(newdata)))
        stop("krls: newdata is missing columns: ",
             paste(setdiff(fi$orig_names, colnames(newdata)),
                   collapse = ", "))
      newdata <- newdata[, fi$orig_names, drop = FALSE]
      # Numeric NA fill with stored training medians.
      if (!is.null(fi$num_medians)) {
        for (jn in names(fi$num_medians)) {
          if (jn %in% names(newdata)) {
            col <- newdata[[jn]]
            if (is.numeric(col) && anyNA(col)) {
              col[is.na(col)] <- fi$num_medians[[jn]]
              newdata[[jn]] <- col
            }
          }
        }
      }
      # Factor NA detection (BUG-010 style).
      fna_report <- character(0)
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
            sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                    jname, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("krls: NA value(s) in factor/character newdata column(s); ",
             "drop the row(s) or impute first.\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      # OOV detection (BUG-004 style).
      oov_report <- list()
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if (is.character(col) || is.factor(col)) {
          col_chr <- as.character(col)
          tr_lev <- fi$xlevels[[jname]]
          bad_rows <- which(!is.na(col_chr) & !(col_chr %in% tr_lev))
          if (length(bad_rows)) {
            bad_lev <- unique(col_chr[bad_rows])
            oov_report[[jname]] <- list(rows = bad_rows, levels = bad_lev)
          }
          newdata[[jname]] <- factor(col, levels = tr_lev)
        }
      }
      if (length(oov_report)) {
        msgs <- vapply(names(oov_report), function(jn) {
          rep <- oov_report[[jn]]
          sprintf("  column %s: %d row(s) (e.g. row %d) with level(s) not seen at fit time: %s",
                  jn, length(rep$rows), rep$rows[1],
                  paste(utils::head(rep$levels, 5), collapse = ", "))
        }, character(1L))
        stop("krls: out-of-vocabulary factor level(s) in newdata.\n",
             paste(msgs, collapse = "\n"), call. = FALSE)
      }
      xnew <- stats::model.matrix(~ ., data = newdata)
      if ("(Intercept)" %in% colnames(xnew))
        xnew <- xnew[, colnames(xnew) != "(Intercept)", drop = FALSE]
      if (!identical(colnames(xnew), fi$expanded_names)) {
        missing_cols <- setdiff(fi$expanded_names, colnames(xnew))
        if (length(missing_cols))
          stop("krls: newdata expansion is missing columns: ",
               paste(missing_cols, collapse = ", "),
               " -- likely an out-of-vocabulary factor level.")
        xnew <- xnew[, fi$expanded_names, drop = FALSE]
      }
      return(as.matrix(xnew))
    }
    if (!is.null(object$terms)) {
      tt <- stats::delete.response(object$terms)
      term_vars <- all.vars(tt)
      fna_report <- character(0)
      for (jn in intersect(term_vars, colnames(newdata))) {
        col <- newdata[[jn]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
            sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                    jn, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("krls: NA value(s) in factor/character newdata column(s).\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      xlev <- if (!is.null(object$xlevels)) object$xlevels else NULL
      mm <- stats::model.matrix(tt, data = newdata, xlev = xlev)
      if ("(Intercept)" %in% colnames(mm))
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      # OOV factor: model.matrix with xlev silently coerces unseen levels
      # to NA in the dummy columns, then na.action=na.omit drops rows.
      # Detect mismatch in row count and stop.
      if (nrow(mm) != nrow(newdata)) {
        stop("krls: newdata expansion dropped rows (likely NA / OOV factor",
             " levels in formula terms). nrow(newdata)=", nrow(newdata),
             ", nrow(mm)=", nrow(mm), ".", call. = FALSE)
      }
      tr_names <- colnames(object$X)
      if (!is.null(tr_names) && !identical(colnames(mm), tr_names)) {
        missing_cols <- setdiff(tr_names, colnames(mm))
        if (length(missing_cols))
          stop("krls: newdata model-matrix expansion is missing columns: ",
               paste(missing_cols, collapse = ", "), call. = FALSE)
        mm <- mm[, tr_names, drop = FALSE]
      }
      return(as.matrix(mm))
    }
    # Plain df, no factor_info, no terms: coerce to matrix.
    return(as.matrix(newdata))
  }
  as.matrix(newdata)
}

## Fix 1: Scale-aware sigma default.
## Returns the median pairwise squared Euclidean distance on the
## standardised predictor matrix Xs (already mean-0, unit-sd per column).
## For n > 500 a 500-row subsample is used for speed (deterministic seed
## 2718 so the result is reproducible within a call).  Floored to
## .Machine$double.eps to prevent degenerate (flat) kernels.
.krls_sigma_anchor <- function(Xs) {
  n <- nrow(Xs)
  if (n > 500L) {
    set.seed(2718L)
    idx <- sample.int(n, 500L)
    Xs_sub <- Xs[idx, , drop = FALSE]
  } else {
    Xs_sub <- Xs
  }
  d2 <- as.numeric(stats::dist(Xs_sub))^2
  max(sqrt(stats::median(d2) * ncol(Xs)), .Machine$double.eps)
}

## Fix 2: Tighter LOO lambda bracket and tolerance.
## tol: changed from 1e-3 * n (n-dependent) to 1e-6 (fixed, 6-digit precision).
## L bracket: changed from linear 0.05-step climb to log-scale x10 climb for
## robustness across a wide range of eigenvalue scales.
.krls_lambdasearch <- function(dvals, V, Vsq, Vty, n_y,
                               L = NULL, U = NULL, tol = NULL,
                               method = "loo", yty = NULL,
                               noisy = FALSE) {
  stopifnot(method %in% c("loo", "gcv"))
  if (method == "gcv") {
    stopifnot(is.numeric(yty), length(yty) == 1L, yty >= 0)
  }
  n <- n_y
  loss_fn <- switch(method,
    loo = function(lam) krls_loo_loss_cpp(dvals, V, Vsq, Vty, lam),
    gcv = function(lam) krls_gcv_loss_cpp(dvals, Vty, yty, n, lam))
  if (is.null(tol)) {
    tol <- 1e-6              ## Fix 2a: was 1e-3 * n
  } else {
    stopifnot(is.numeric(tol), length(tol) == 1, tol > 0)
  }
  if (is.null(U)) {
    U <- n
    while (sum(dvals / (dvals + U)) < 1) U <- U - 1
  } else {
    stopifnot(is.numeric(U), length(U) == 1, U > 0)
  }
  if (is.null(L)) {
    q  <- which.min(abs(dvals - (max(dvals) / 1000)))
    L  <- .Machine$double.eps
    while (sum(dvals / (dvals + L)) > q) L <- L * 10  ## Fix 2b: was L + 0.05
  } else {
    stopifnot(is.numeric(L), length(L) == 1, L >= 0)
  }
  gr <- 0.381966
  X1 <- L + gr * (U - L)
  X2 <- U - gr * (U - L)
  S1 <- loss_fn(X1)
  S2 <- loss_fn(X2)
  if (noisy) {
    cat("L:", L, "X1:", X1, "X2:", X2, "U:", U, "S1:", S1, "S2:", S2, "\n")
  }
  while (abs(S1 - S2) > tol) {
    if (S1 < S2) {
      U  <- X2; X2 <- X1
      X1 <- L + gr * (U - L)
      S2 <- S1
      S1 <- loss_fn(X1)
    } else {
      L  <- X1; X1 <- X2
      X2 <- U - gr * (U - L)
      S1 <- S2
      S2 <- loss_fn(X2)
    }
    if (noisy) {
      cat("L:", L, "X1:", X1, "X2:", X2, "U:", U,
          "S1:", S1, "S2:", S2, "\n")
    }
  }
  if (S1 < S2) X1 else X2
}

.krls_fd_binary <- function(object) {
  d <- ncol(object$X); n <- nrow(object$X)
  lu <- function(x) length(unique(x))
  binidx <- which(apply(object$X, 2L, lu) == 2L)
  if (length(binidx) == 0L) return(object)
  est  <- se <- matrix(NA_real_, nrow = 1L, ncol = length(binidx))
  diffs_store <- matrix(NA_real_, nrow = n, ncol = length(binidx))
  for (i in seq_along(binidx)) {
    X1 <- X0 <- object$X
    X1[, binidx[i]] <- max(X1[, binidx[i]])
    X0[, binidx[i]] <- min(X0[, binidx[i]])
    Xall <- rbind(X1, X0)
    h <- matrix(rep(c(1 / n, -(1 / n)), each = n), ncol = 1L)
    pout <- predict.krls_rr(object, newdata = Xall, se.fit = TRUE)
    est[1L, i] <- as.numeric(t(h) %*% pout$fit)
    se[1L, i]  <- as.numeric(sqrt(t(h) %*% pout$vcov.fit %*% h)) * sqrt(2)
    diffs_store[, i] <- pout$fit[1:n] - pout$fit[(n + 1):(2 * n)]
  }
  object$derivatives[, binidx]        <- diffs_store
  object$avgderivatives[, binidx]     <- est
  object$var.avgderivatives[, binidx] <- se^2
  object$binaryindicator[, binidx]    <- TRUE
  object
}

#' @export
print.krls_rr <- function(x, ...) {
  cat("Kernel Regularized Least Squares (KRLS)\n")
  cat("  n =", nrow(x$X), "  p =", ncol(x$X), "\n")
  cat("  sigma =", signif(x$sigma, 4),
      "  lambda =", signif(x$lambda, 4),
      "  R^2 =", signif(x$R2, 4), "\n")
  if (!is.null(x$weights)) {
    cat("  weighted: yes (mean(w)=1 normalisation)\n")
  }
  if (!is.null(x$autotune)) {
    cat("  Autotune: sigma=", signif(x$autotune$winner, 4),
        " (grid: ", paste(signif(x$autotune$grid, 3), collapse = ", "),
        ")\n", sep = "")
  }
  if (!is.null(x$cv)) {
    cat("  CV: ", x$cv$nfold, "-fold x ", x$cv$ncross,
        " cross  (grid: ", length(x$cv$lambda.grid), " lambda)",
        if (isTRUE(x$cv$cv.1se)) "  [1-SE rule]" else "",
        "\n", sep = "")
  }
  if (!is.null(x$boot)) {
    cat("  Bagging: n.boot = ", x$boot$n.boot, " replicate(s)\n",
        sep = "")
  }
  if (!is.null(x$avgderivatives)) {
    cat("\nAverage Marginal Effects:\n")
    print(setNames(as.vector(x$avgderivatives),
                   colnames(x$avgderivatives)))
  }
  invisible(x)
}

#' @export
summary.krls_rr <- function(object, ...) {
  out <- list()
  out$call_info <- c(n = nrow(object$X), p = ncol(object$X),
                     sigma = object$sigma, lambda = object$lambda,
                     R2 = object$R2)
  out$weighted <- !is.null(object$weights)
  out$autotune <- object$autotune
  out$cv <- object$cv
  out$boot_n <- if (!is.null(object$boot)) object$boot$n.boot else NULL
  if (!is.null(object$avgderivatives)) {
    se <- sqrt(as.numeric(object$var.avgderivatives))
    avg <- as.numeric(object$avgderivatives)
    tval <- avg / se
    pval <- 2 * pt(-abs(tval), df = nrow(object$X) - 1L)
    coefs <- cbind(Estimate   = avg,
                   `Std. Err` = se,
                   `t value`  = tval,
                   `Pr(>|t|)` = pval)
    rownames(coefs) <- colnames(object$avgderivatives)
    out$avg_eff <- coefs
    quart <- apply(object$derivatives, 2L, quantile,
                   probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    out$quartiles <- quart
  }
  class(out) <- "summary.krls_rr"
  out
}

#' @export
print.summary.krls_rr <- function(x, ...) {
  cat("Kernel Regularized Least Squares (KRLS)\n")
  ci <- x$call_info
  cat(sprintf("  n = %d  p = %d  sigma = %s  lambda = %s  R^2 = %s\n",
              as.integer(ci["n"]), as.integer(ci["p"]),
              signif(ci["sigma"], 4),
              signif(ci["lambda"], 4),
              signif(ci["R2"], 4)))
  if (isTRUE(x$weighted)) {
    cat("  weighted: yes (mean(w)=1 normalisation)\n")
  }
  if (!is.null(x$autotune)) {
    cat("  Autotune: sigma=", signif(x$autotune$winner, 4), "\n",
        sep = "")
  }
  if (!is.null(x$cv)) {
    cat("  CV: ", x$cv$nfold, "-fold x ", x$cv$ncross,
        " cross\n", sep = "")
  }
  if (!is.null(x$boot_n)) {
    cat("  Bagging: n.boot = ", x$boot_n, "\n", sep = "")
  }
  if (!is.null(x$avg_eff)) {
    cat("\nAverage Marginal Effects:\n")
    printCoefmat(x$avg_eff, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nQuartiles of Pointwise Marginal Effects:\n")
    print(x$quartiles)
  }
  invisible(x)
}

#' Diagnostic plot for KRLS fits
#'
#' Four-panel diagnostic display (residuals vs fitted, Q-Q on
#' standardised residuals, scale-location, residuals vs leverage with
#' Cook contours) in a 2x2 grid by default. Panels 4 (Cook's distance)
#' and 6 (Cook vs leverage) are available via `which = 1:6`. Modelled
#' on `stats::plot.lm()` and `plot.ares()`.
#'
#' Leverage uses the KRLS hat matrix
#' \deqn{H = K (K + \lambda I)^{-1} = V \mathrm{diag}(d/(d+\lambda)) V'}
#' so `h_i = sum_k V[i,k]^2 * d_k / (d_k + lambda)`. The effective
#' degrees of freedom is `sum(d/(d+lambda))`. Residual dispersion is
#' `sqrt(RSS / (n - effdf))` (or `out$varmod$sigma_hat` when present).
#'
#' @param x A fitted `krls_rr` object.
#' @param which Which panels to draw; integer subset of `1:6`. Defaults
#'   to `c(1L, 2L, 3L, 5L)`.
#' @param caption Panel titles (length 6).
#' @param sub.caption Optional sub-caption shown across the figure.
#' @param main Per-panel title.
#' @param ask Activate interactive panel paging.
#' @param id.n Number of extreme points to label per panel.
#' @param labels.id Optional vector of point labels.
#' @param cex.id Label cex.
#' @param qqline Draw the Q-Q line.
#' @param cook.levels Cook's-distance contour levels.
#' @param add.smooth Add LOESS smoother overlay.
#' @param label.pos Label position codes.
#' @param panel Per-panel display function (defaults to LOESS smoother).
#' @param ... Passed to underlying `plot` calls.
#' @return Invisibly returns `x`.
#' @export
plot.krls_rr <- function(x,
                         which = c(1L, 2L, 3L, 5L),
                         caption = list("Residuals vs Fitted", "Normal Q-Q",
                                        "Scale-Location",
                                        "Cook's distance",
                                        "Residuals vs Leverage",
                                        expression("Cook's dist vs Leverage  " *
                                                   h[ii] / (1 - h[ii]))),
                         panel = if (add.smooth) graphics::panel.smooth
                                 else graphics::points,
                         sub.caption = NULL, main = "",
                         ask = prod(graphics::par("mfcol")) <
                                 length(which) &&
                               grDevices::dev.interactive(),
                         ...,
                         id.n = 3L,
                         labels.id = NULL,
                         cex.id = 0.75,
                         qqline = TRUE,
                         cook.levels = c(0.5, 1.0),
                         add.smooth = getOption("add.smooth", TRUE),
                         label.pos = c(4, 2)) {
  if (!inherits(x, "krls"))
    stop("plot.krls_rr: 'x' must be a 'krls' object.")
  which <- as.integer(which)
  if (any(!which %in% 1:6))
    stop("plot.krls_rr: `which` must be a subset of 1:6.")
  show <- rep(FALSE, 6L); show[which] <- TRUE

  yhat <- as.numeric(x$fitted)
  yv   <- as.numeric(x$y)
  rraw <- yv - yhat
  n <- length(yhat)
  if (n == 0L) stop("plot.krls_rr: fit has zero observations.")
  w <- if (!is.null(x$weights)) as.numeric(x$weights) else rep(1, n)

  # Hat matrix from K's eigen-decomp. We re-decompose K here (cheap
  # for n < 2000; large-n users can pre-compute via the fit's stored
  # K if they want to skip plot).
  K <- x$K
  if (is.null(K))
    stop("plot.krls_rr: fit has no stored K matrix; cannot compute leverages.")
  eo <- krls_eig_cpp(K)
  dvals <- as.numeric(eo$values)
  V     <- eo$vectors
  lam   <- x$lambda
  shrink <- dvals / (dvals + lam)
  # Hat diag: h_i = sum_k V[i,k]^2 * shrink_k
  h <- as.numeric((V * V) %*% shrink)
  h <- pmin(pmax(h, 0), 1 - .Machine$double.eps)

  effdf <- sum(shrink)
  rdf <- max(n - effdf, 1)
  sigma_hat <- if (!is.null(x$varmod) && !is.null(x$varmod$sigma_hat))
    x$varmod$sigma_hat
  else sqrt(sum(w * rraw^2) / rdf)

  rstd <- rraw * sqrt(w) / (sigma_hat * sqrt(1 - h))
  rstd[!is.finite(rstd)] <- NA_real_

  cook <- (rstd^2 / max(effdf, 1)) * h / (1 - h)
  cook[!is.finite(cook)] <- NA_real_

  if (is.null(labels.id)) labels.id <- as.character(seq_len(n))
  extrm <- function(v, k = id.n) {
    if (k < 1L) return(integer(0))
    finite <- which(is.finite(v))
    if (!length(finite)) return(integer(0))
    ord <- finite[order(-abs(v[finite]))]
    ord[seq_len(min(k, length(ord)))]
  }

  one_fig <- all(graphics::par("mfcol") == c(1L, 1L))
  if (one_fig && length(which) > 1L) {
    op <- graphics::par(mfrow = c(2L, 2L), oma = c(0, 0, 2, 0),
                        mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
  } else if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op), add = TRUE)
  }

  if (is.null(sub.caption)) {
    cl <- if (!is.null(x$call))
      paste(deparse(x$call, width.cutoff = 75L), collapse = " ")
    else "krls fit"
    sub.caption <- if (nchar(cl) > 90L)
      paste0(substring(cl, 1L, 87L), "...") else cl
  }

  if (show[1L]) {
    ylim <- range(rraw, na.rm = TRUE)
    if (id.n > 0L) ylim <- grDevices::extendrange(r = ylim, f = 0.08)
    graphics::plot(yhat, rraw, xlab = "Fitted values",
                   ylab = "Residuals",
                   main = main, ylim = ylim, type = "n", ...)
    panel(yhat, rraw, ...)
    graphics::abline(h = 0, lty = 3, col = "gray")
    idx <- extrm(rraw)
    if (length(idx))
      graphics::text(yhat[idx], rraw[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[1L]], side = 3, line = 0.25, cex = 0.85)
  }
  if (show[2L]) {
    qq <- stats::qqnorm(rstd, main = main,
                        ylab = "Standardized residuals", ...)
    if (qqline) stats::qqline(rstd, lty = 3, col = "gray")
    idx <- extrm(rstd)
    if (length(idx))
      graphics::text(qq$x[idx], qq$y[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[2L]], side = 3, line = 0.25, cex = 0.85)
  }
  if (show[3L]) {
    sqrtabs <- sqrt(abs(rstd))
    graphics::plot(yhat, sqrtabs, xlab = "Fitted values",
                   ylab = expression(sqrt(abs(`Standardized residuals`))),
                   main = main, type = "n", ...)
    panel(yhat, sqrtabs, ...)
    idx <- extrm(sqrtabs)
    if (length(idx))
      graphics::text(yhat[idx], sqrtabs[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[3L]], side = 3, line = 0.25, cex = 0.85)
  }
  if (show[4L]) {
    ylim <- c(0, max(cook, na.rm = TRUE) * 1.075)
    graphics::plot(seq_len(n), cook, type = "h", main = main,
                   xlab = "Obs. number", ylab = "Cook's distance",
                   ylim = ylim, ...)
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(idx, cook[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[4L]], side = 3, line = 0.25, cex = 0.85)
  }
  if (show[5L]) {
    xlim <- c(0, max(h, na.rm = TRUE) * 1.05)
    ylim <- range(rstd, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-3, 3)
    graphics::plot(h, rstd, xlim = xlim, ylim = ylim, xlab = "Leverage",
                   ylab = "Standardized residuals", main = main,
                   type = "n", ...)
    panel(h, rstd, ...)
    graphics::abline(h = 0, v = 0, lty = 3, col = "gray")
    if (length(cook.levels)) {
      hh <- seq.int(0.001, max(h, na.rm = TRUE), length.out = 101L)
      hh <- hh[hh < 1 & hh > 0]
      for (cl_lvl in cook.levels) {
        rr <- sqrt(cl_lvl * max(effdf, 1) * (1 - hh) / hh)
        graphics::lines(hh, rr, lty = 2, col = "red")
        graphics::lines(hh, -rr, lty = 2, col = "red")
      }
    }
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(h[idx], rstd[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[5L]], side = 3, line = 0.25, cex = 0.85)
  }
  if (show[6L]) {
    hsc <- h / (1 - h)
    graphics::plot(hsc, cook,
                   xlab = expression(h[ii] / (1 - h[ii])),
                   ylab = "Cook's distance", main = main, type = "n", ...)
    panel(hsc, cook, ...)
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(hsc[idx], cook[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[6L]], side = 3, line = 0.25, cex = 0.85)
  }

  if (!is.null(sub.caption) && length(which) > 1L && one_fig)
    graphics::mtext(sub.caption, outer = TRUE, cex = 0.9, line = 0.5)

  invisible(x)
}
