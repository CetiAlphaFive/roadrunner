# R/plda.R

#' Penalized linear discriminant analysis
#'
#' Penalized Fisher's linear discriminant (Witten & Tibshirani 2011) with L1 or
#' fused-lasso penalties, multi-class support, and built-in cross-validation.
#'
#' @param x Numeric predictor matrix, or a formula for the formula interface.
#' @param ... Passed to methods.
#' @return An object of S3 class `"plda"`.
#' @export
plda <- function(x, ...) UseMethod("plda")

#' @param y Factor (or coercible) class label of length `nrow(x)`.
#' @param K Number of discriminant vectors (`<= G-1`). If supplied, `K` is
#'   fixed; if `NULL` (default), `K` is chosen by cross-validation along with
#'   `lambda`. When `autotune = FALSE`, `NULL` means `G-1`.
#' @param lambda Magnitude penalty. Required when `autotune = FALSE`.
#' @param penalty `"L1"` (default) or `"fused"`.
#' @param lambda2 Fused-lasso difference penalty (used when `penalty = "fused"`).
#' @param autotune If `TRUE` (default), cross-validate `lambda` (and `K`).
#' @param nfold CV folds.
#' @param lambda_grid Optional CV grid.
#' @param maxit MM solver maximum iterations.
#' @param tol MM solver convergence tolerance.
#' @param nthreads Integer. Number of worker threads for the cross-validation
#'   autotune (`autotune = TRUE`). `0` (the default, taken from
#'   `getOption("roadrunner.nthreads")` when set) uses
#'   `RcppParallel::defaultNumThreads()`. The CV fold loop is parallelised with
#'   fixed per-fold output slots and a serial-order reduction, so the fitted
#'   discriminants are byte-identical regardless of `nthreads`. Ignored when
#'   `autotune = FALSE` (a single fit is already sequential).
#' @rdname plda
#' @export
plda.default <- function(x, y, K = NULL, lambda = NULL,
                         penalty = c("L1", "fused"), lambda2 = NULL,
                         autotune = TRUE, nfold = 5L, lambda_grid = NULL,
                         maxit = 100L, tol = 1e-6,
                         nthreads = 0L, ...) {
  # Save and restore the global RNG state so plda() does not leave .Random.seed
  # behind when the caller had none — every Rcpp export (plda_wcsd_cpp,
  # plda_fit_cpp, plda_project_cpp) creates .Random.seed via RNGScope.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
  }
  penalty <- match.arg(penalty)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("plda: `x` must be numeric.", call. = FALSE)
  if (anyNA(x)) stop("plda: `x` contains missing values (NA). Remove or impute before fitting.", call. = FALSE)
  if (ncol(x) < 1L) stop("plda: `x` must have at least one column.", call. = FALSE)
  y <- droplevels(as.factor(y))
  if (anyNA(y))
    stop("plda: `y` contains missing class labels (NA). Remove rows with missing labels before fitting.", call. = FALSE)
  if (nrow(x) != length(y)) stop("plda: nrow(x) must equal length(y).", call. = FALSE)
  classes <- levels(y)
  G <- length(classes)
  if (G < 2L) stop("plda: need at least two classes.", call. = FALSE)
  yint <- as.integer(y)
  # LDA needs within-class variance, so every class must have >= 2 members.
  # A singleton class also makes stratified CV folds impossible.
  class_counts <- tabulate(yint, nbins = G)
  if (any(class_counts < 2L)) {
    bad <- classes[class_counts < 2L]
    stop(sprintf("plda: every class needs at least 2 observations; too few in: %s.",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  # An explicit K is honoured under autotune (only lambda is tuned); K = NULL
  # lets cross-validation choose K alongside lambda.
  tune_K <- is.null(K)
  if (is.null(K)) K <- G - 1L
  if (K < 1L || K > G - 1L)
    stop(sprintf("plda: K must be in 1..%d (G-1).", G - 1L), call. = FALSE)
  lam2 <- if (is.null(lambda2)) 0 else lambda2
  pen_code <- if (penalty == "fused") 1L else 0L

  # Validate nfold
  nfold <- as.integer(nfold)
  if (is.na(nfold) || nfold < 2L)
    stop("plda: `nfold` must be an integer >= 2.", call. = FALSE)
  if (nfold > nrow(x))
    stop(sprintf("plda: `nfold` (%d) cannot exceed nrow(x) (%d).", nfold, nrow(x)),
         call. = FALSE)
  # BUG-021 (v0.0.0.9055): the stratified fold builder places each class
  # into folds 1..length(class) via sample(rep_len(seq_len(nfold), ...)),
  # which means when min(class_count) < nfold some folds receive zero
  # test observations for that class. Those folds short-circuit with
  # err_slot = 0, biasing the CV error toward 0 and selecting the wrong
  # lambda. Reject upfront with a clear message recommending the
  # largest valid nfold.
  if (autotune) {
    min_class <- min(class_counts)
    if (nfold > min_class)
      stop(sprintf(
        paste0("plda: `nfold` (%d) exceeds the smallest class count ",
               "(%d). Stratified CV cannot place every class into ",
               "every fold; set nfold = %d (or smaller)."),
        nfold, min_class, min_class),
        call. = FALSE)
  }

  # Resolve nthreads. `0` (or anything <= 0) means "use the package default".
  # Mirror krls.default: bare 0L in the signature; option lookup + resolution
  # in the body so `getOption("roadrunner.nthreads")` is honoured at call time.
  # Determinism invariant: the parallel CV harness (plda_cv_inner_cpp) writes
  # each fold into a private slot and reduces in fixed fold order, so the
  # fitted discriminants are byte-identical for any nthreads.
  nthreads <- as.integer(nthreads)
  if (length(nthreads) != 1L || is.na(nthreads) || nthreads < 0L)
    stop("plda: `nthreads` must be a single non-negative integer.", call. = FALSE)
  if (nthreads == 0L) {
    opt_nt <- getOption("roadrunner.nthreads", 0L)
    nthreads <- if (!is.null(opt_nt) && opt_nt > 0L) as.integer(opt_nt) else 0L
  }
  nthreads_eff <- if (nthreads <= 0L)
    RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  # Validate lambda2 / warn if irrelevant
  if (penalty == "L1" && !is.null(lambda2) && lambda2 != 0) {
    warning("plda: `lambda2` is ignored when `penalty = 'L1'`.", call. = FALSE)
  }
  if (penalty == "fused") {
    if (!is.numeric(lam2) || length(lam2) != 1L || !is.finite(lam2) || lam2 < 0)
      stop("plda: `lambda2` must be a single non-negative finite number when `penalty = 'fused'`.", call. = FALSE)
  }

  nthreads_used <- 1L
  if (autotune) {
    cv <- .plda_cv(x, yint, G, K, penalty, pen_code, lam2, nfold,
                   lambda_grid, maxit, tol, nthreads_eff, tune_K)
    lambda <- cv$lambda; K <- cv$K
    nthreads_used <- cv$nthreads_used
  } else if (is.null(lambda)) {
    stop("plda: supply `lambda` or use autotune = TRUE.", call. = FALSE)
  } else {
    cv <- NULL
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop("plda: `lambda` must be a single positive finite number.", call. = FALSE)

  eng <- plda_fit_cpp(x, yint, G, K, lambda, lam2, pen_code,
                      as.integer(maxit), tol)
  structure(list(discrim = eng$discrim, mu = eng$mu, sdw = eng$sdw,
                 cmeans = eng$cmeans, classes = classes,
                 K = K, lambda = lambda, lambda2 = lam2, penalty = penalty,
                 cv = cv, nthreads_used = nthreads_used, call = match.call()),
            class = "plda")
}

#' @keywords internal
#' @noRd
# Default data-driven log-spaced lambda grid.
# `hi` is the largest standardized between-class mean across features.
.plda_lambda_grid <- function(x, yint, G, n = 12L) {
  S <- plda_wcsd_cpp(x, yint, G)
  S[S < 1e-12] <- 1
  # standardize: centre by global mean, scale by within-class sd
  xs <- scale(x, center = TRUE, scale = S)
  # per-class means in standardized space (rows = classes)
  cm <- t(vapply(seq_len(G),
                 function(g) colMeans(xs[yint == g, , drop = FALSE]),
                 numeric(ncol(x))))
  # upper lambda: largest standardized between-class signal across features
  hi <- max(abs(cm)) + 1e-8
  exp(seq(log(hi * 1e-3), log(hi), length.out = n))
}

#' @keywords internal
#' @noRd
# k-fold CV over the lambda grid; picks the lambda (and, when tune_K = TRUE,
# the K in 1..K) with lowest mean misclassification error. When tune_K = FALSE
# the supplied K is held fixed and only lambda is tuned. Saves and restores the
# global RNG state so the caller's random-number stream is unaffected (same
# pattern as .ares_autotune and .krls_lambdasearch in this package).
#
# Parallelism: the fold assignment is computed here in R (set.seed(0) +
# sample) and passed into plda_cv_inner_cpp, which runs the fold loop with TBB
# — one fold per task, each fold writing its misclassification counts into a
# private slot, reduced in fixed fold order. The result is therefore byte-
# identical for any `nthreads` (the roadrunner determinism invariant). `nthreads`
# is a trailing argument with a default of 1L so existing positional callers
# (e.g. tests calling .plda_cv directly) remain valid and deterministic.
.plda_cv <- function(x, yint, G, K, penalty, pen_code, lam2, nfold,
                     lambda_grid, maxit, tol, nthreads = 1L, tune_K = TRUE) {
  # Save and restore the global RNG state so this call does not leak a seed
  # change to the caller — consistent with .ares_autotune and .krls_lambdasearch.
  # MUST be the very first statements: .plda_lambda_grid calls plda_wcsd_cpp
  # (an Rcpp export) which creates .Random.seed via RNGScope/GetRNGstate(), so
  # the existence check must run before any C++ is touched.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
  }
  if (is.null(lambda_grid)) lambda_grid <- .plda_lambda_grid(x, yint, G)
  n <- nrow(x)
  set.seed(0L)
  # Class-stratified fold assignment: each class is spread evenly across folds
  # so every training fold (the nfold-1 complement) retains every class. With
  # the >= 2-per-class guarantee from plda.default this means no training fold
  # can lose a class. Still deterministic — set.seed(0L) fixes the sampling.
  folds <- integer(n)
  for (g in seq_len(G)) {
    ig <- which(yint == g)
    folds[ig] <- sample(rep_len(seq_len(nfold), length(ig)))
  }

  res <- plda_cv_inner_cpp(x, as.integer(yint), as.integer(folds),
                           as.integer(nfold), G, K, as.numeric(lambda_grid),
                           lam2, pen_code, as.integer(maxit), tol,
                           as.integer(nthreads))
  err <- res$err / n
  if (!isTRUE(res$ok))
    warning("plda: minorize-maximize criterion decreased in one or more CV folds; ",
            "results may be unreliable. Consider a smaller lambda grid or larger maxit.",
            call. = FALSE)
  # When the user supplied K explicitly (tune_K = FALSE), tune lambda only:
  # restrict the (nlam x K) error matrix to the single K-th column so the K
  # choice is fixed at the supplied value.
  if (!tune_K) err <- err[, K, drop = FALSE]
  best <- which(err == min(err), arr.ind = TRUE)[1, ]
  best_K <- if (tune_K) as.integer(best[2]) else K
  list(lambda = lambda_grid[best[1]], K = best_K,
       grid = lambda_grid, errors = err[, best[2]],
       nthreads_used = as.integer(res$nthreads_used))
}

#' @param formula A model formula; response on the left, predictors on the right.
#' @param data A data frame (or environment) containing the formula variables.
#' @rdname plda
#' @export
plda.formula <- function(formula, data = NULL, ...) {
  if (is.null(data)) data <- environment(formula)
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  x <- stats::model.matrix(attr(mf, "terms"), mf)
  icpt <- match("(Intercept)", colnames(x), nomatch = 0L)
  if (icpt > 0L) x <- x[, -icpt, drop = FALSE]
  fit <- plda.default(x, y, ...)
  fit$call <- match.call()
  fit$terms <- attr(mf, "terms")
  fit$xlevels <- stats::.getXlevels(attr(mf, "terms"), mf)
  fit
}
