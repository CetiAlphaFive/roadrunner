# logreg -- binary logistic regression fitter
#
# File layout:
#   1. logreg()           -- S3 generic
#   2. logreg.formula()   -- formula -> matrix dispatch (mirrors ols.formula)
#   3. logreg.default()   -- main entry; validation, C++ IRLS engine call,
#                            optional bagging.
#
# Conventions:
#   - The C++ engine is reached via two Rcpp exports: logreg_fit_cpp and
#     logreg_vcov_cpp (see src/logistic.cpp).
#   - An intercept column is included by default (formula path via
#     model.matrix; default path prepends a column of ones).
#   - The fit minimises the binomial deviance by iteratively reweighted
#     least squares (Fisher scoring). Non-convergence -- the classic
#     symptom of (quasi-)perfect separation -- is reported with a warning
#     rather than an error.
#   - Determinism contract: bagging uses a serial bootstrap loop, so fits
#     are byte-identical across thread counts. Any change that breaks this
#     invariant is a bug.

# ============================================================================
#  S3 generic + formula method
# ============================================================================

#' Binary logistic regression
#'
#' Fits a binary logistic regression model by iteratively reweighted least
#' squares (Fisher scoring). The IRLS step is a weighted economy-QR solve
#' in C++ (Armadillo/LAPACK), mirroring the [ols()] engine. The response
#' must be binary; it may be supplied as a 0/1 numeric vector, a logical
#' vector, or a two-level factor.
#'
#' Two interfaces are provided. The formula method takes a model formula
#' and a data frame. The default method takes a numeric predictor matrix
#' `x` and a binary response `y`.
#'
#' Optional bagging (`n.boot > 0`) refits the model on bootstrap-weight
#' resamples; `predict()` then returns the bag-mean probability. Bagging
#' is a serial loop and is byte-identical across thread counts.
#'
#' When the IRLS iteration fails to converge within `maxit` steps -- the
#' usual signature of (quasi-)complete separation -- the fit is returned
#' with `converged = FALSE` and a warning is emitted.
#'
#' @param x A numeric matrix or data frame of predictors, or a model
#'   formula when calling the formula method.
#' @param y A binary response of length `nrow(x)`: a 0/1 numeric vector, a
#'   logical vector, or a factor with exactly two levels (the second
#'   level is treated as the "success" class).
#' @param data A data frame. Used only by the formula method.
#' @param subset Optional row-subsetting vector for the formula method
#'   (integer indices or logical), passed through to
#'   [stats::model.frame()]. Ignored by the default method.
#' @param na.action Strategy for missing values, passed to
#'   [stats::model.frame()] for the formula method. Default
#'   [stats::na.omit()]. For the default method, rows with any `NA` in
#'   `x` or `y` are dropped with a warning.
#' @param weights Optional non-negative numeric vector of length
#'   `nrow(x)` of prior weights for weighted logistic regression. Default
#'   `NULL` (unweighted). Negative, `NA`, `NaN`, or `Inf` weights are
#'   rejected; zero weights are rejected (drop the row instead).
#' @param intercept Logical; include an intercept term. Default `TRUE`.
#'   Honoured by the default method; the formula method follows the
#'   formula's own intercept specification.
#' @param n.boot Number of bootstrap replicate fits for bagging.
#'   Default `0` (no bagging). When `> 0`, `predict()` averages the
#'   fitted probabilities over the replicates plus the central fit.
#' @param seed Optional integer seed for the bagging bootstrap draws.
#'   Pass an integer for reproducible bagging; `NULL` (default) uses the
#'   current RNG state. The user's RNG stream is restored on exit when a
#'   seed is supplied.
#' @param maxit Maximum number of IRLS iterations. Default `25`.
#' @param tol Convergence tolerance on the relative change in deviance
#'   between IRLS iterations. Default `1e-8`.
#' @param ... Currently ignored.
#' @return An object of class `"logreg"`: a list with the fitted
#'   `coefficients`, `fitted.values` (probabilities), `linear.predictors`,
#'   deviance and working `residuals`, `deviance`, `null.deviance`, `aic`,
#'   `df.residual`, `df.null`, the classical variance-covariance matrix
#'   `vcov`, the IRLS `iter` count and `converged` flag, the numerical
#'   `rank`, `hatvalues`, `(X'WX)^-1` as `XtWXinv`, the design matrix `X`
#'   and the (numeric 0/1) response `y`, echoed `weights`, the `call`,
#'   and (when requested) a `$boot` field.
#' @references
#' McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*,
#' 2nd ed. Chapman and Hall.
#'
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' *Journal of Econometrics* 29(3):305-325.
#' @examples
#' df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
#' fit <- logreg(y ~ wt + hp, data = df)
#' print(fit)
#' predict(fit, df[1:3, ], type = "response")
#' @export
logreg <- function(x, ...) UseMethod("logreg")

#' @rdname logreg
#' @export
logreg.formula <- function(x, data = NULL, subset = NULL,
                           na.action = stats::na.omit, weights = NULL, ...,
                           y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()

  # Reject offset() terms: logreg has no offset concept; silently
  # absorbing one into the linear predictor would be wrong.
  pre_tt <- stats::terms(formula, data = data)
  if (length(attr(pre_tt, "offset"))) {
    stop("logreg: offset() terms are not supported in the formula. ",
         "Subtract the offset from the linear predictor before calling ",
         "logreg() if needed.", call. = FALSE)
  }

  # Build the model frame. Weights pass through model.frame so the
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
    stop("logreg: response variable is missing from formula/data.")
  w <- stats::model.weights(mf)
  mm <- stats::model.matrix(formula, mf)

  out <- logreg.default(x = mm, y = yv, weights = w,
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

#' @rdname logreg
#' @export
logreg.default <- function(x, y, weights = NULL, intercept = TRUE,
                           n.boot = 0L, seed = NULL,
                           maxit = 25L, tol = 1e-8, ...) {
  cl <- match.call()

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
    stop("logreg: x must be a numeric matrix or data frame.")
  if (length(y) != nrow(x))
    stop("logreg: length(y) (", length(y), ") must equal nrow(x) (",
         nrow(x), ").")
  storage.mode(x) <- "double"

  # ---- Binary response coercion -------------------------------------------
  # Accept a 0/1 numeric vector, a logical vector, or a two-level factor.
  # The non-NA values must reduce to exactly two distinct classes; the
  # response is stored as numeric 0/1 with the "second" class as 1.
  y_levels <- NULL
  if (is.factor(y)) {
    lv <- levels(y)
    if (length(lv) != 2L)
      stop("logreg: a factor response must have exactly two levels; got ",
           length(lv), ".")
    y_levels <- lv
    y <- as.numeric(y) - 1
  } else if (is.logical(y)) {
    y_levels <- c("FALSE", "TRUE")
    y <- as.numeric(y)
  } else if (is.numeric(y)) {
    uy <- sort(unique(y[!is.na(y)]))
    if (length(uy) > 2L || (length(uy) == 2L && !all(uy == c(0, 1))) ||
        (length(uy) == 1L && !uy %in% c(0, 1)))
      stop("logreg: a numeric response must be coded 0/1; observed ",
           "value(s): ", paste(utils::head(uy, 5L), collapse = ", "),
           if (length(uy) > 5L) ", ..." else "", ".")
    y <- as.numeric(y)
  } else {
    stop("logreg: y must be a 0/1 numeric vector, a logical vector, or a ",
         "two-level factor.")
  }

  # ---- Observation weights (validation) -----------------------------------
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("logreg: weights must be numeric.")
    if (length(weights) != length(y))
      stop("logreg: length(weights) (", length(weights),
           ") must equal length(y) (", length(y), ").")
    if (any(is.na(weights)) || any(!is.finite(weights[!is.na(weights)])))
      stop("logreg: weights contain NA / NaN / Inf; all weights must be ",
           "finite.")
  }

  # ---- Missing-value handling (default method) ----------------------------
  # The default method drops any row with NA in x, y, or weights. The
  # formula method already applies its own na.action upstream.
  na_mask <- !stats::complete.cases(x) | is.na(y)
  if (!is.null(weights)) na_mask <- na_mask | is.na(weights)
  if (any(na_mask)) {
    warning("logreg: dropped ", sum(na_mask),
            " row(s) with missing values.", call. = FALSE)
    keep <- !na_mask
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(weights)) weights <- weights[keep]
  }
  if (any(!is.finite(x)))
    stop("logreg: x contains NaN / Inf values; aborting fit.")

  # Re-check the response is genuinely binary after NA removal: a column
  # that collapses to a single class carries no logistic information.
  uy <- unique(y)
  if (length(uy) < 2L)
    stop("logreg: response has only one class after dropping missing ",
         "values; logistic regression needs both classes present.")

  # Reject non-positive weights now that NA rows are gone. Zero weights
  # would bias the residual df, so we require strictly positive weights.
  if (!is.null(weights) && any(weights <= 0))
    stop("logreg: weights must be strictly positive. To omit a row, drop",
         " it from x/y; weights = 0 is not a substitute.")

  n <- nrow(x)
  if (n < 1L) stop("logreg: no observations left to fit.")

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
  if (p < 1L) stop("logreg: x must have at least one column.")
  if (n <= p)
    stop("logreg: need more observations than columns (n = ", n,
         ", p = ", p, ").")
  has_int <- isTRUE(intercept) || "(Intercept)" %in% colnames(x)

  # ---- IRLS controls ------------------------------------------------------
  maxit <- as.integer(maxit)
  if (is.na(maxit) || maxit < 1L)
    stop("logreg: maxit must be a positive integer.")
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0)
    stop("logreg: tol must be a single positive number.")

  # ---- Weight vector for the engine ---------------------------------------
  # The engine expects an explicit prior-weight vector; unweighted passes
  # ones.
  w_vec <- if (is.null(weights)) rep(1, n) else as.numeric(weights)

  # ---- C++ IRLS fit -------------------------------------------------------
  # BUG-020 (v0.0.0.9055): pass has_intercept so the null-deviance
  # baseline uses mu0 = 0.5 for no-intercept formulas (y ~ 0 + x), where
  # the canonical null model is eta = 0 -> mu = 0.5. Without this, the
  # engine always used the weighted mean of y and disagreed with
  # stats::glm() while df.null was already adjusted on the R side.
  eng <- logreg_fit_cpp(x, y, w_vec, maxit, tol, has_int)
  coef <- as.numeric(eng$coefficients)
  names(coef) <- colnames(x)

  if (!isTRUE(eng$converged))
    warning("logreg: IRLS did not converge in ", maxit, " iteration(s). ",
            "This often indicates (quasi-)perfect separation; ",
            "coefficient estimates may be unreliable.", call. = FALSE)

  # ---- Residuals ----------------------------------------------------------
  # Deviance residuals: sign(y - mu) * sqrt(unit deviance). Working
  # residuals: (y - mu) / (mu (1 - mu)), the IRLS working response minus
  # the linear predictor.
  mu <- as.numeric(eng$fitted)
  eta <- as.numeric(eng$eta)
  dev_unit <- numeric(n)
  for (i in seq_len(n)) {
    yi <- y[i]; mi <- mu[i]
    term <- 0
    if (yi > 0) term <- term + yi * log(yi / mi)
    if (yi < 1) term <- term + (1 - yi) * log((1 - yi) / (1 - mi))
    dev_unit[i] <- 2 * w_vec[i] * term
  }
  dev_resid <- sign(y - mu) * sqrt(pmax(dev_unit, 0))
  work_resid <- (y - mu) / (mu * (1 - mu))

  # AIC for the binomial model: deviance + 2 * (number of parameters).
  # With binary y the saturated log-likelihood is 0, so the model
  # log-likelihood is -deviance / 2 and aic = deviance + 2 p.
  aic <- eng$deviance + 2 * p

  out <- list(
    coefficients      = coef,
    fitted.values     = mu,
    linear.predictors = eta,
    residuals         = dev_resid,
    working.residuals = work_resid,
    working.weights   = as.numeric(eng$working),
    deviance          = eng$deviance,
    null.deviance     = eng$`null.deviance`,
    aic               = aic,
    df.residual       = n - p,
    df.null           = n - as.integer(has_int),
    rank              = p,
    iter              = eng$iter,
    converged         = isTRUE(eng$converged),
    hatvalues         = as.numeric(eng$hatdiag),
    XtWXinv           = eng$XtWXinv,
    vcov              = eng$XtWXinv,
    X                 = x,
    y                 = y,
    weights           = weights,
    y.levels          = y_levels,
    intercept         = has_int,
    n                 = n,
    p                 = p,
    maxit             = maxit,
    tol               = tol,
    call              = cl,
    factor_info       = factor_info
  )
  dimnames(out$vcov) <- list(colnames(x), colnames(x))
  dimnames(out$XtWXinv) <- list(colnames(x), colnames(x))

  # ---- Bagging ------------------------------------------------------------
  # Bootstrap-weight resamples: each replicate draws row multiplicities
  # ~ Multinomial(n, 1/n), multiplies them into the prior weights, and
  # refits. Serial loop -> trivially deterministic across thread counts.
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
      # A resample collapsing to one class is skipped (NA column).
      yb <- y[pos]
      if (length(unique(yb)) < 2L) next
      f_b <- tryCatch(
        logreg_fit_cpp(x[pos, , drop = FALSE], yb, w_b[pos], maxit, tol,
                       has_int),
        error = function(e) NULL)
      if (!is.null(f_b))
        boot_coef[, b] <- as.numeric(f_b$coefficients)
    }
    rownames(boot_coef) <- colnames(x)
    out$boot <- list(coefficients = boot_coef, n.boot = n_boot)
  }

  class(out) <- "logreg"
  out
}
