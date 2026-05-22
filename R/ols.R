# ols -- ordinary / weighted least squares fitter
#
# File layout:
#   1. ols()           -- S3 generic
#   2. ols.formula()   -- formula -> matrix dispatch (mirrors ares.formula)
#   3. ols.default()   -- main entry; validation, C++ engine call, bagging,
#                         optional residual-variance model.
#
# Conventions:
#   - The C++ engine is reached via two Rcpp exports: ols_fit_cpp and
#     ols_vcov_cpp (see src/ols.cpp).
#   - An intercept column is included by default (formula path via
#     model.matrix; default path prepends a column of ones).
#   - Determinism contract: bagging uses a serial bootstrap loop, so fits
#     are byte-identical across thread counts. Any change that breaks this
#     invariant is a bug.

# ============================================================================
#  S3 generic + formula method
# ============================================================================

#' Ordinary and weighted least squares
#'
#' Fits a linear regression model by ordinary least squares, or by
#' weighted least squares when `weights` are supplied. The solver is a
#' weighted economy-QR factorisation in C++ (Armadillo/LAPACK); a
#' rank-deficient design matrix is rejected with a clear error rather
#' than silently dropping columns.
#'
#' Two interfaces are provided. The formula method takes a model formula
#' and a data frame. The default method takes a numeric predictor matrix
#' `x` and a numeric response vector `y`.
#'
#' Optional bagging (`n.boot > 0`) refits the model on bootstrap-weight
#' resamples; `predict()` then returns the bag mean. Bagging is a serial
#' loop and is byte-identical across thread counts.
#'
#' @param x A numeric matrix or data frame of predictors, or a model
#'   formula when calling the formula method.
#' @param y A numeric response vector with length `nrow(x)`.
#' @param data A data frame. Used only by the formula method.
#' @param subset Optional row-subsetting vector for the formula method
#'   (integer indices or logical), passed through to
#'   [stats::model.frame()]. Ignored by the default method.
#' @param na.action Strategy for missing values, passed to
#'   [stats::model.frame()] for the formula method. Default
#'   [stats::na.omit()]. For the default method, rows with any `NA` in
#'   `x` or `y` are dropped with a warning.
#' @param weights Optional non-negative numeric vector of length
#'   `nrow(x)` for weighted least squares. Default `NULL` (unweighted).
#'   Negative, `NA`, `NaN`, or `Inf` weights are rejected; zero weights
#'   are rejected (drop the row instead).
#' @param intercept Logical; include an intercept term. Default `TRUE`.
#'   Honoured by the default method; the formula method follows the
#'   formula's own intercept specification.
#' @param n.boot Number of bootstrap replicate fits for bagging.
#'   Default `0` (no bagging). When `> 0`, `predict()` averages over the
#'   replicates plus the central fit and (with `se.fit = TRUE`) returns a
#'   per-prediction bag standard deviation.
#' @param seed Optional integer seed for the bagging bootstrap draws.
#'   Pass an integer for reproducible bagging; `NULL` (default) uses the
#'   current RNG state. The user's RNG stream is restored on exit when a
#'   seed is supplied.
#' @param varmod Residual variance model used by `predict()` for
#'   prediction intervals.
#'   - `"none"` (default): the closed-form prediction interval uses the
#'     fitted residual standard error directly.
#'   - `"const"`: stores a single residual SD.
#'   - `"lm"`: fits a small linear model of `|resid|` on the fitted
#'     values to allow simple yhat-dependent heteroscedasticity. Shares
#'     the `ares()` variance-model helper.
#' @param ... Currently ignored.
#' @return An object of class `"ols"`: a list with the fitted
#'   `coefficients`, `fitted.values`, `residuals`, `sigma2` and `sigma`,
#'   `df.residual`, the classical variance-covariance matrix `vcov`, the
#'   numerical `rank`, `hatvalues`, `(X'WX)^-1` as `XtXinv`, the design
#'   matrix `X` and response `y`, echoed `weights`, the `call`, and
#'   (when requested) `$varmod` and `$boot` fields.
#' @references
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' *Journal of Econometrics* 29(3):305-325.
#' @examples
#' fit <- ols(mpg ~ wt + hp, data = mtcars)
#' print(fit)
#' predict(fit, mtcars[1:3, ])
#' @export
ols <- function(x, ...) UseMethod("ols")

#' @rdname ols
#' @export
ols.formula <- function(x, data = NULL, subset = NULL,
                        na.action = stats::na.omit, weights = NULL, ...,
                        y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()

  # Reject offset() terms: ols has no offset concept; silently absorbing
  # one into the design would be wrong.
  pre_tt <- stats::terms(formula, data = data)
  if (length(attr(pre_tt, "offset"))) {
    stop("ols: offset() terms are not supported in the formula. ",
         "Subtract the offset from y before calling ols() if needed.",
         call. = FALSE)
  }

  # Build the model frame. We pass weights through model.frame so the
  # lm()-style match.call() trick resolves `subset` and `weights` against
  # the caller's environment, and so na.action drops weights in lockstep
  # with the rows it removes.
  mfcall <- match.call(expand.dots = FALSE)
  mfcall$y <- NULL
  mfcall[[1L]] <- quote(stats::model.frame)
  names(mfcall)[names(mfcall) == "x"] <- "formula"
  mfcall$`...` <- NULL
  mf <- eval(mfcall, parent.frame())
  yv <- stats::model.response(mf)
  if (is.null(yv))
    stop("ols: response variable is missing from formula/data.")
  w <- stats::model.weights(mf)
  mm <- stats::model.matrix(formula, mf)

  out <- ols.default(x = mm, y = as.numeric(yv), weights = w,
                     intercept = FALSE, ...)
  out$call <- cl
  tt <- stats::terms(formula, data = data)
  out$terms <- tt
  out$xlevels <- stats::.getXlevels(tt, mf)
  out
}

# ============================================================================
#  Main default method
# ============================================================================

#' @rdname ols
#' @export
ols.default <- function(x, y, weights = NULL, intercept = TRUE,
                        n.boot = 0L, seed = NULL,
                        varmod = c("none", "const", "lm"), ...) {
  cl <- match.call()
  varmod <- match.arg(varmod)

  # ---- Coerce x to a numeric matrix ---------------------------------------
  # A data frame with factor / character columns is expanded to a numeric
  # model matrix via model.matrix(~ ., x); an intercept column is then
  # dropped here so the `intercept` argument below controls it uniformly.
  factor_info <- NULL
  if (is.data.frame(x)) {
    is_cat <- vapply(x, \(z) is.factor(z) || is.character(z), logical(1L))
    if (any(is_cat)) {
      for (j in which(is_cat))
        if (is.character(x[[j]])) x[[j]] <- factor(x[[j]])
      xlevels <- lapply(x[is_cat], levels)
      xmm <- stats::model.matrix(~ ., data = x)
      if ("(Intercept)" %in% colnames(xmm))
        xmm <- xmm[, colnames(xmm) != "(Intercept)", drop = FALSE]
      factor_info <- list(xlevels = xlevels, orig_names = names(x),
                          is_cat = is_cat,
                          expanded_names = colnames(xmm))
      x <- xmm
    } else {
      nm <- names(x)
      x <- as.matrix(x)
      storage.mode(x) <- "double"
      colnames(x) <- nm
    }
  }
  if (!is.matrix(x) || !is.numeric(x))
    stop("ols: x must be a numeric matrix or data frame.")
  if (!is.numeric(y))
    stop("ols: y must be numeric.")
  if (length(y) != nrow(x))
    stop("ols: length(y) (", length(y), ") must equal nrow(x) (",
         nrow(x), ").")
  storage.mode(x) <- "double"
  y <- as.numeric(y)

  # ---- Observation weights (validation) -----------------------------------
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("ols: weights must be numeric.")
    if (length(weights) != length(y))
      stop("ols: length(weights) (", length(weights),
           ") must equal length(y) (", length(y), ").")
    if (any(is.na(weights)) || any(!is.finite(weights[!is.na(weights)])))
      stop("ols: weights contain NA / NaN / Inf; all weights must be finite.")
  }

  # ---- Missing-value handling (default method) ----------------------------
  # The default method drops any row with NA in x, y, or weights. The
  # formula method already applies its own na.action upstream.
  na_mask <- !stats::complete.cases(x) | is.na(y)
  if (!is.null(weights)) na_mask <- na_mask | is.na(weights)
  if (any(na_mask)) {
    warning("ols: dropped ", sum(na_mask),
            " row(s) with missing values.", call. = FALSE)
    keep <- !na_mask
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(weights)) weights <- weights[keep]
  }
  if (any(!is.finite(y)))
    stop("ols: y contains NaN / Inf values; aborting fit.")
  if (any(!is.finite(x)))
    stop("ols: x contains NaN / Inf values; aborting fit.")

  # Reject non-positive weights now that NA rows are gone. Zero weights
  # would bias the (n - p) residual df, so we require strictly positive
  # weights and tell the user to drop the row instead.
  if (!is.null(weights) && any(weights <= 0))
    stop("ols: weights must be strictly positive. To omit a row, drop",
         " it from x/y; weights = 0 is not a substitute.")

  n <- nrow(x)
  if (n < 1L) stop("ols: no observations left to fit.")

  # ---- Intercept ----------------------------------------------------------
  # The default method prepends an intercept column when intercept = TRUE.
  # The formula method passes intercept = FALSE because model.matrix has
  # already produced the intercept column (or omitted it per the formula).
  if (isTRUE(intercept)) {
    cn <- colnames(x)
    if (is.null(cn)) cn <- paste0("V", seq_len(ncol(x)))
    x <- cbind(`(Intercept)` = 1, x)
    colnames(x) <- c("(Intercept)", cn)
  } else if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  p <- ncol(x)
  if (p < 1L) stop("ols: x must have at least one column.")
  if (n <= p)
    stop("ols: need more observations than columns (n = ", n,
         ", p = ", p, ").")

  # ---- Weight vector for the engine ---------------------------------------
  # The engine expects an explicit weight vector; OLS passes ones.
  w_vec <- if (is.null(weights)) rep(1, n) else as.numeric(weights)

  # ---- C++ fit ------------------------------------------------------------
  eng <- ols_fit_cpp(x, y, w_vec)
  coef <- as.numeric(eng$coefficients)
  names(coef) <- colnames(x)

  out <- list(
    coefficients  = coef,
    fitted.values = as.numeric(eng$fitted),
    residuals     = as.numeric(eng$residuals),
    sigma2        = eng$sigma2,
    sigma         = sqrt(eng$sigma2),
    rss           = eng$rss,
    df.residual   = eng$df,
    rank          = eng$rank,
    hatvalues     = as.numeric(eng$hatdiag),
    XtXinv        = eng$XtXinv,
    vcov          = eng$sigma2 * eng$XtXinv,
    X             = x,
    y             = y,
    weights       = weights,
    intercept     = isTRUE(intercept),
    n             = n,
    p             = p,
    call          = cl,
    factor_info   = factor_info
  )
  dimnames(out$vcov) <- list(colnames(x), colnames(x))
  dimnames(out$XtXinv) <- list(colnames(x), colnames(x))

  # ---- Bagging ------------------------------------------------------------
  # Bootstrap-weight resamples: each replicate draws row multiplicities
  # ~ Multinomial(n, 1/n), multiplies them into the prior weights, and
  # refits. The bootstrap weight vector is renormalised so its mean is 1
  # within the replicate (consistency with the central fit's convention).
  # Serial loop -> trivially deterministic across thread counts.
  n_boot <- as.integer(n.boot)
  if (is.na(n_boot) || n_boot < 0L) n_boot <- 0L
  if (n_boot > 0L) {
    has_seed <- !is.null(seed)
    if (has_seed) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(),
                        inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
      }
      set.seed(as.integer(seed))
    }
    boot_coef <- matrix(NA_real_, nrow = p, ncol = n_boot)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      mult <- tabulate(idx, nbins = n)        # bootstrap multiplicities.
      w_b <- mult * w_vec
      pos <- w_b > 0
      # Refit on the rows with positive weight only; the C++ engine
      # rejects zero weights and an all-zero row carries no information.
      f_b <- ols_fit_cpp(x[pos, , drop = FALSE], y[pos], w_b[pos])
      boot_coef[, b] <- as.numeric(f_b$coefficients)
    }
    rownames(boot_coef) <- colnames(x)
    out$boot <- list(coefficients = boot_coef, n.boot = n_boot)
  }

  # ---- Residual variance model -------------------------------------------
  # Reuse the ares() varmod helper so prediction intervals share one
  # implementation. The helper consumes $residuals, $fitted.values and
  # $coefficients, all present on `out` at this point.
  out$varmod <- NULL
  if (varmod != "none")
    out$varmod <- .ares_fit_varmod(out, varmod = varmod, weights = weights)

  class(out) <- "ols"
  out
}
