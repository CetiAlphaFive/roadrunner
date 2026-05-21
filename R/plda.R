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
#' @param K Number of discriminant vectors (`<= G-1`). Default `G-1`.
#' @param lambda Magnitude penalty. Required when `autotune = FALSE`.
#' @param penalty `"L1"` (default) or `"fused"`.
#' @param lambda2 Fused-lasso difference penalty (used when `penalty = "fused"`).
#' @param autotune If `TRUE` (default), cross-validate `lambda` (and `K`).
#' @param nfold CV folds.
#' @param lambda_grid Optional CV grid.
#' @param maxit MM solver maximum iterations.
#' @param tol MM solver convergence tolerance.
#' @rdname plda
#' @export
plda.default <- function(x, y, K = NULL, lambda = NULL,
                         penalty = c("L1", "fused"), lambda2 = NULL,
                         autotune = TRUE, nfold = 5L, lambda_grid = NULL,
                         maxit = 100L, tol = 1e-6, ...) {
  penalty <- match.arg(penalty)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("plda: `x` must be numeric.", call. = FALSE)
  if (anyNA(x)) stop("plda: `x` contains missing values (NA). Remove or impute before fitting.", call. = FALSE)
  if (ncol(x) < 1L) stop("plda: `x` must have at least one column.", call. = FALSE)
  y <- as.factor(y)
  if (anyNA(y))
    stop("plda: `y` contains missing class labels (NA). Remove rows with missing labels before fitting.", call. = FALSE)
  if (nrow(x) != length(y)) stop("plda: nrow(x) must equal length(y).", call. = FALSE)
  classes <- levels(y)
  G <- length(classes)
  if (G < 2L) stop("plda: need at least two classes.", call. = FALSE)
  yint <- as.integer(y)
  if (is.null(K)) K <- G - 1L
  if (K < 1L || K > G - 1L)
    stop(sprintf("plda: K must be in 1..%d (G-1).", G - 1L), call. = FALSE)
  lam2 <- if (is.null(lambda2)) 0 else lambda2
  pen_code <- if (penalty == "fused") 1L else 0L

  # Validate lambda2 / warn if irrelevant
  if (penalty == "L1" && !is.null(lambda2) && lambda2 != 0) {
    warning("plda: `lambda2` is ignored when `penalty = 'L1'`.", call. = FALSE)
  }
  if (penalty == "fused") {
    if (!is.numeric(lam2) || length(lam2) != 1L || !is.finite(lam2) || lam2 < 0)
      stop("plda: `lambda2` must be a single non-negative finite number when `penalty = 'fused'`.", call. = FALSE)
  }

  if (autotune) {
    cv <- .plda_cv(x, yint, G, K, penalty, pen_code, lam2, nfold,
                   lambda_grid, maxit, tol)
    lambda <- cv$lambda; K <- cv$K
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
                 cmeans = eng$cmeans, cw = eng$cw, classes = classes,
                 K = K, lambda = lambda, lambda2 = lam2, penalty = penalty,
                 cv = cv, call = match.call()),
            class = "plda")
}

# @keywords internal
# Default data-driven log-spaced lambda grid.
.plda_lambda_grid <- function(x, yint, G, n = 12L) {
  S <- plda_wcsd_cpp(x, yint, G)
  S[S < 1e-12] <- 1
  xs <- scale(x, center = TRUE, scale = FALSE) %*% diag(1 / S)
  hi <- max(abs(colMeans(xs))) + 1e-8
  exp(seq(log(hi * 1e-3), log(hi), length.out = n))
}

# @keywords internal
# k-fold CV over (lambda grid) x (1..K); picks the lambda/K with lowest mean
# misclassification error. Saves and restores the global RNG state so the
# caller's random-number stream is unaffected (same pattern as .ares_autotune
# and .krls_lambdasearch in this package).
.plda_cv <- function(x, yint, G, K, penalty, pen_code, lam2, nfold,
                     lambda_grid, maxit, tol) {
  if (is.null(lambda_grid)) lambda_grid <- .plda_lambda_grid(x, yint, G)
  n <- nrow(x)
  # Save and restore the global RNG state so this call does not leak a seed
  # change to the caller — consistent with .ares_autotune and .krls_lambdasearch.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(0L)
  folds <- sample(rep_len(seq_len(nfold), n))
  err <- matrix(0, length(lambda_grid), K)
  for (f in seq_len(nfold)) {
    tr <- folds != f; te <- !tr
    for (li in seq_along(lambda_grid)) {
      eng <- plda_fit_cpp(x[tr, , drop = FALSE], yint[tr], G, K,
                          lambda_grid[li], lam2, pen_code,
                          as.integer(maxit), tol)
      sc  <- plda_project_cpp(x[te, , drop = FALSE], eng$mu, eng$sdw, eng$discrim)
      csc <- eng$cmeans %*% eng$discrim
      for (k in seq_len(K)) {
        d2 <- outer(rowSums(sc[, 1:k, drop = FALSE]^2),
                    rowSums(csc[, 1:k, drop = FALSE]^2), `+`) -
              2 * sc[, 1:k, drop = FALSE] %*% t(csc[, 1:k, drop = FALSE])
        pred <- max.col(-d2, ties.method = "first")
        err[li, k] <- err[li, k] + sum(pred != yint[te])
      }
    }
  }
  err <- err / n
  best <- which(err == min(err), arr.ind = TRUE)[1, ]
  list(lambda = lambda_grid[best[1]], K = as.integer(best[2]),
       grid = lambda_grid, errors = err[, best[2]])
}

#' @param formula A model formula; response on the left, predictors on the right.
#' @param data A data frame (or environment) containing the formula variables.
#' @rdname plda
#' @export
plda.formula <- function(formula, data = NULL, ...) {
  if (is.null(data)) data <- environment(formula)
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(attr(mf, "terms"), mf)
  icpt <- match("(Intercept)", colnames(x), nomatch = 0L)
  if (icpt > 0L) x <- x[, -icpt, drop = FALSE]
  fit <- plda.default(x, y, ...)
  fit$call <- match.call()
  fit$terms <- attr(mf, "terms")
  fit$xlevels <- .getXlevels(attr(mf, "terms"), mf)
  fit
}
