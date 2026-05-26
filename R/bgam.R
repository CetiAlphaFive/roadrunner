# bgam -- component-wise P-spline gradient boosting
#
# File layout:
#   1. bgam()             -- S3 generic
#   2. bgam.formula()     -- formula -> matrix dispatch (mirrors ols.formula)
#   3. bgam.default()     -- main entry: validation, basis construction,
#                            lambda calibration, boosting loop, optional
#                            autotune CV and bagging.
#   4. .bgam_lambda_from_df()  -- calibrate lambda_j from df_target
#   5. .bgam_cv_mstop()        -- k-fold CV to select mstop_opt
#   6. .bgam_bagged_fit()      -- n.boot bootstrap-weight resamples
#
# Conventions:
#   - All internal helpers prefixed with `.bgam_`.
#   - C++ engines: bgam_bspline_basis_cpp, bgam_diff_penalty_cpp,
#     bgam_prefactor_cpp, bgam_boost_cpp, bgam_predict_cpp
#     (see src/bgam.cpp).
#   - Determinism contract: TBB parallel_reduce over predictors resolves
#     ties by lowest j index; nthreads=1 vs nthreads=N is byte-identical.
#     Bagging uses a serial RNG loop -- byte-identical across thread counts.

# ============================================================================
#  S3 generic + formula method
# ============================================================================

#' Component-wise P-spline gradient boosting
#'
#' Fits a smooth, strictly additive model via component-wise functional
#' gradient boosting with penalised B-spline (P-spline) base-learners, one
#' per predictor column. At each of `mstop` boosting iterations the algorithm
#' selects the single base-learner that most reduces the current loss (squared
#' error for gaussian, binomial deviance for binomial), updates only that
#' component by a shrinkage-scaled Cholesky-backed ridge solve, and repeats.
#' The result is a smooth GAM with implicit variable selection: predictors
#' that are never selected retain a zero coefficient vector, and
#' `$selection_frequency` (fraction of iterations each predictor was chosen)
#' serves as a natural importance metric.
#'
#' Each base-learner is a cubic (or general degree-\code{d}) B-spline with
#' \code{nknots} interior knots placed at equally-spaced quantiles of the
#' predictor, penalised by a \code{dpen}-th order difference-penalty matrix
#' \eqn{D^T D}{D'D} scaled by a ridge parameter \eqn{\lambda_j}{lambda_j}
#' calibrated so the effective degrees of freedom of the base-learner equals
#' \code{df_target}.  Near-constant predictors fall back automatically to
#' unpenalised linear base-learners (with a warning).
#'
#' `bgam()` complements the rest of the roadrunner lineup: it is more
#' interpretable than `krls()` (explicit additive components, no kernel
#' matrix), handles non-linearity without manual basis specification unlike
#' `ols()`, and targets smooth monotonic / additive signals where `ares()`'s
#' piecewise-linear hinge functions may be less efficient.  The complexity per
#' fit is O(n * p * K) per iteration (K ~ `nknots + degree + 1`), making it
#' practical for large-n, moderately large-p settings where `krls()` becomes
#' memory-bound.
#'
#' The formula method builds a standard model frame (factors expanded via
#' `model.matrix`, formula's `-1` disables the intercept initialisation),
#' stores `$terms` and `$xlevels` for `predict()`, and dispatches to the
#' default method.  The default method takes a numeric matrix directly.
#'
#' @param x A numeric matrix or data frame of predictors (default method), or
#'   a model formula (formula method).
#' @param y Response vector of length `nrow(x)`.  For `family = "binomial"`,
#'   must be binary: 0/1 numeric, `logical`, or a two-level `factor` (second
#'   level = 1).
#' @param data A data frame (formula method only).
#' @param subset Optional integer or logical row-subsetting vector passed to
#'   [stats::model.frame()] (formula method only).
#' @param na.action NA-handling function passed to [stats::model.frame()].
#'   Default [stats::na.omit()].  The default method drops rows with any `NA`
#'   in `x` and emits a warning; rows with `NA` in `y` are always rejected.
#' @param weights Optional strictly positive numeric weight vector of length
#'   `nrow(x)`.  For `"gaussian"`, weights scale the pseudo-residual sum of
#'   squares used for base-learner selection and the Cholesky prefactor.  For
#'   `"binomial"`, weights scale the IRLS working weights
#'   \eqn{w_i^{(\text{eff})} = w_i \mu_i (1 - \mu_i)}{w_i * mu_i * (1 - mu_i)}.
#'   Default `NULL` (unweighted).
#' @param family Response family: `"gaussian"` (default, squared-error loss)
#'   or `"binomial"` (logistic / binomial deviance, binary response).
#' @param nu Boosting shrinkage / learning rate.  Must be in `(0, 1]`.
#'   Smaller values require more iterations but produce smoother solution
#'   paths.  Default `0.1`.
#' @param mstop Number of boosting iterations to run.  If `NULL` and
#'   `autotune = FALSE`, defaults to `mstop_max`.  If `autotune = TRUE`,
#'   `mstop` is set to the CV-selected `mstop_opt`.
#' @param autotune If `TRUE` (default), selects `mstop` by k-fold
#'   cross-validation over the grid `seq(10, mstop_max, by = 10)` using MSE
#'   (gaussian) or binomial deviance (binomial) as the criterion.  A
#'   successive-halving early-stop rule terminates the grid scan once CV loss
#'   rises more than 2% of its range above the running best.
#' @param mstop_max Maximum iterations for the CV grid (and the default when
#'   `autotune = FALSE` and `mstop = NULL`).  Reduced automatically to
#'   `max(50, 2 * n)` when `n < 100`.  Default `300L`.
#' @param nfold Number of CV folds.  Default `5L`.
#' @param seed.cv Integer seed for CV fold construction.  Default `NULL`.
#'   Binomial CV folds are stratified by response class.
#' @param nknots Number of interior B-spline knots per predictor.  Knots are
#'   placed at equally-spaced quantiles of the predictor column.  Default
#'   `20L`, giving `K_j = nknots + degree + 1 = 24` columns per base-learner
#'   at defaults.
#' @param degree B-spline polynomial degree (1 to 5).  Default `3L` (cubic
#'   B-splines).
#' @param dpen Order of the difference penalty (1 to `degree`).  Default `2L`
#'   (second-order penalty, analogous to a thin-plate-spline roughness penalty;
#'   Eilers & Marx 1996).
#' @param df_target Target effective degrees of freedom per base-learner,
#'   i.e. `tr(S_j(lambda_j))`.  Must be `> 1`.  Default `4`.  The ridge
#'   penalty `lambda_j` is calibrated per predictor via `stats::uniroot`.
#' @param lambda_method `"df"` (default): calibrate `lambda_j` from
#'   `df_target` per predictor.  `"fixed"`: supply `lambda_j` directly via
#'   `lambda_fixed`.
#' @param lambda_fixed Numeric scalar or length-p vector of lambda values.
#'   Required when `lambda_method = "fixed"`; ignored otherwise.
#' @param unpenalized Optional character vector of predictor column names to
#'   enter as linear (unpenalised, no-spline) base-learners.  Useful for
#'   known-linear covariates or binary dummies.
#' @param intercept If `TRUE` (default), initialise the additive predictor at
#'   `mean(y)` (gaussian) or `logit(clamp(mean(y), 0.001, 0.999))` (binomial).
#'   If `FALSE`, initialise at 0.
#' @param n.boot Number of bootstrap replicates for bagging.  Default `0L`
#'   (no bagging).  When `> 0`, `predict()` returns the bag mean across the
#'   central fit and surviving replicates.  Hyperparameters are frozen at
#'   the central-fit values; no re-autotune per replicate.
#' @param seed Optional integer RNG seed for the bagging loop.  Default
#'   `NULL`.
#' @param nthreads TBB thread count for the parallel predictor scan per
#'   boosting iteration.  `0` = use all available TBB threads.  `1L`
#'   (default) = serial execution.  Results are byte-identical across all
#'   `nthreads` values.
#' @param ... Passed to the default method from the formula method; otherwise
#'   currently unused.
#'
#' @return An S3 object of class `"bgam"`, a named list with the following
#'   key fields (full contract: 25+ fields):
#'   \describe{
#'     \item{`coefficients`}{List of p accumulated (nu-scaled) beta vectors,
#'       one per predictor.  Zero-vector if the predictor was never selected.}
#'     \item{`fitted.values`}{Fitted values on the response scale
#'       (probabilities for binomial; identical to `linear.predictors` for
#'       gaussian).  For bagged fits, this is the bag mean.}
#'     \item{`linear.predictors`}{Fitted values on the link scale.}
#'     \item{`residuals`}{Residuals on the response scale (gaussian:
#'       `y - fitted`; binomial: deviance residuals).}
#'     \item{`selection_path`}{Integer vector of length `mstop`, 1-based
#'       predictor index selected at each iteration.}
#'     \item{`selection_frequency`}{Named numeric vector of length p; fraction
#'       of iterations each predictor was selected.  Sums to 1.}
#'     \item{`loss_path`}{Numeric vector of length `mstop`; training loss
#'       (MSE or deviance) after each iteration.}
#'     \item{`mstop`}{Integer; iterations run (equals `mstop_opt` when
#'       `autotune = TRUE`).}
#'     \item{`cv`}{List with `$grid`, `$cv_loss`, `$mstop_opt` when
#'       `autotune = TRUE`; `NULL` otherwise.}
#'     \item{`boot`}{List with `$n.boot` and `$fits` when `n.boot > 0`;
#'       `NULL` otherwise.}
#'     \item{`base_learners`}{List of p per-predictor lists, each containing
#'       `$knots`, `$K`, `$degree`, `$dpen`, `$lambda`, `$range`, `$chol`,
#'       `$B_train`.}
#'   }
#'   See also: `predict.bgam`, `print.bgam`, `summary.bgam`, `plot.bgam`.
#'
#' @references
#' Buehlmann, P. and Yu, B. (2003). Boosting with the L2 loss: Regression and
#' classification. *Journal of the American Statistical Association*,
#' 98(462):324--339. \doi{10.1198/016214503000125}
#'
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11(2):89--121.
#' \doi{10.1214/ss/1038425655}
#'
#' Hofner, B., Mayr, A., Robinzonov, N. and Schmid, M. (2014). Model-based
#' boosting in R: A hands-on tutorial using the R package mboost.
#' *Computational Statistics*, 29(1--2):3--35.
#' \doi{10.1007/s00180-012-0382-5}
#'
#' Schmid, M. and Hothorn, T. (2008). Boosting additive models using
#' component-wise P-splines. *Computational Statistics and Data Analysis*,
#' 53(2):298--311. \doi{10.1016/j.csda.2008.09.009}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 4), n, 4)
#' y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
#'
#' # Fit with fixed mstop (fast, no CV)
#' fit <- bgam(x, y, mstop = 100, autotune = FALSE, nthreads = 1)
#' print(fit)
#' predict(fit, x[1:3, ])
#'
#' # Formula interface
#' df <- data.frame(y = y, x)
#' fit2 <- bgam(y ~ ., data = df, mstop = 50, autotune = FALSE)
#' head(fit2$selection_frequency)
#' @export
bgam <- function(x, ...) UseMethod("bgam")

#' @rdname bgam
#' @export
bgam.formula <- function(x, data = NULL, subset = NULL,
                         na.action = stats::na.omit, weights = NULL, ...,
                         y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()

  # Reject offset() terms
  pre_tt <- stats::terms(formula, data = data)
  if (length(attr(pre_tt, "offset")))
    stop("bgam: offset() terms are not supported in the formula.",
         call. = FALSE)

  # Check for - 1 (no-intercept formula)
  has_intercept_formula <- attr(pre_tt, "intercept") == 1L

  # Build model frame
  mfcall <- match.call(expand.dots = FALSE)
  mfcall$y <- NULL
  mfcall[[1L]] <- quote(stats::model.frame)
  names(mfcall)[names(mfcall) == "x"] <- "formula"
  mfcall$`...` <- NULL
  mf <- eval(mfcall, parent.frame())
  yv <- stats::model.response(mf)
  if (is.null(yv))
    stop("bgam: response variable is missing from formula/data.")
  w <- stats::model.weights(mf)
  mm <- stats::model.matrix(formula, mf)

  # Drop intercept column from mm; bgam.default handles initialization
  if ("(Intercept)" %in% colnames(mm))
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  out <- bgam.default(x = mm, y = yv, weights = w,
                      intercept = has_intercept_formula, ...)
  out$call <- cl
  tt <- stats::terms(formula, data = data)
  out$terms <- tt
  out$xlevels <- stats::.getXlevels(tt, mf)
  out
}

# ============================================================================
#  Main default method
# ============================================================================

#' @rdname bgam
#' @export
bgam.default <- function(x, y,
                         family        = c("gaussian", "binomial"),
                         nu            = 0.1,
                         mstop         = NULL,
                         autotune      = TRUE,
                         mstop_max     = 300L,
                         nfold         = 5L,
                         seed.cv       = NULL,
                         nknots        = 20L,
                         degree        = 3L,
                         dpen          = 2L,
                         df_target     = 4,
                         lambda_method = c("df", "fixed"),
                         lambda_fixed  = NULL,
                         unpenalized   = NULL,
                         intercept     = TRUE,
                         weights       = NULL,
                         n.boot        = 0L,
                         seed          = NULL,
                         nthreads      = 1L,
                         ...) {
  cl <- match.call()
  # Validate family before match.arg so the error message is user-friendly.
  if (!is.null(family) && length(family) == 1L &&
      !family %in% c("gaussian", "binomial")) {
    stop(sprintf("bgam: family must be one of 'gaussian', 'binomial'; got '%s'.",
                 family), call. = FALSE)
  }
  family        <- match.arg(family)
  lambda_method <- match.arg(lambda_method)

  # ---- Coerce x to numeric matrix -----------------------------------------
  if (is.data.frame(x)) {
    nm_orig <- names(x)
    is_cat <- vapply(x, function(z) is.factor(z) || is.character(z), logical(1L))
    if (any(is_cat)) {
      for (j in which(is_cat))
        if (is.character(x[[j]])) x[[j]] <- factor(x[[j]])
      xmm <- stats::model.matrix(~ ., data = x)
      if ("(Intercept)" %in% colnames(xmm))
        xmm <- xmm[, colnames(xmm) != "(Intercept)", drop = FALSE]
      x <- xmm
    } else {
      x <- as.matrix(x)
      colnames(x) <- nm_orig
    }
  }
  if (!is.matrix(x) || !is.numeric(x))
    stop("bgam: x must be a numeric matrix or data frame.")
  if (ncol(x) == 0L)
    stop("bgam: x must have at least 1 column.", call. = FALSE)
  storage.mode(x) <- "double"

  # ---- Binary response coercion (binomial only) ---------------------------
  if (identical(family, "binomial")) {
    if (is.factor(y)) {
      lv <- levels(y)
      if (length(lv) != 2L)
        stop("bgam: y must be binary (0/1, logical, or 2-level factor) for ",
             "family='binomial'.", call. = FALSE)
      y <- as.numeric(y) - 1
    } else if (is.logical(y)) {
      y <- as.numeric(y)
    } else if (is.numeric(y)) {
      uy <- sort(unique(y[!is.na(y)]))
      if (!((length(uy) <= 2L) &&
            (length(uy) == 0L ||
             (length(uy) == 1L && uy %in% c(0, 1)) ||
             (length(uy) == 2L && all(uy == c(0, 1))))))
        stop("bgam: y must be binary (0/1, logical, or 2-level factor) for ",
             "family='binomial'.", call. = FALSE)
      y <- as.numeric(y)
    } else {
      stop("bgam: y must be binary (0/1, logical, or 2-level factor) for ",
           "family='binomial'.", call. = FALSE)
    }
  } else {
    if (!is.numeric(y))
      stop("bgam: y must be numeric for family='gaussian'.", call. = FALSE)
    y <- as.numeric(y)
  }

  # ---- Basic validation ---------------------------------------------------
  if (length(y) != nrow(x))
    stop("bgam: length(y) must equal nrow(x).", call. = FALSE)
  if (anyNA(y))
    stop("bgam: NA values in y are not allowed.", call. = FALSE)
  if (any(!is.finite(y)))
    stop("bgam: Non-finite values in y are not allowed.", call. = FALSE)

  # ---- NA handling in x (drop rows with warning) --------------------------
  na_mask <- !stats::complete.cases(x)
  na_idx <- NULL
  if (any(na_mask)) {
    warning("bgam: dropped ", sum(na_mask),
            " row(s) with NA in x.", call. = FALSE)
    na_idx <- which(na_mask)
    keep <- !na_mask
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(weights)) weights <- weights[keep]
  }

  n <- nrow(x)
  p_raw <- ncol(x)

  if (n < 10L)
    stop("bgam: need at least 10 observations.", call. = FALSE)

  # ---- Weights validation -------------------------------------------------
  if (!is.null(weights)) {
    if (length(weights) != n)
      stop("bgam: length(weights) must equal nrow(x).", call. = FALSE)
    if (any(weights <= 0 | is.na(weights) | !is.finite(weights)))
      stop("bgam: weights must be strictly positive and finite.", call. = FALSE)
    weights <- as.numeric(weights)
  }

  # ---- Parameter validation -----------------------------------------------
  nknots    <- as.integer(nknots)
  degree    <- as.integer(degree)
  dpen      <- as.integer(dpen)
  mstop_max <- as.integer(mstop_max)
  nfold     <- as.integer(nfold)
  nthreads  <- as.integer(nthreads)

  if (is.na(nknots) || nknots < 1L)
    stop("bgam: nknots must be >= 1.", call. = FALSE)
  if (is.na(degree) || degree < 1L || degree > 5L)
    stop("bgam: degree must be between 1 and 5.", call. = FALSE)
  if (is.na(dpen) || dpen < 1L || dpen >= degree + 1L)
    stop("bgam: dpen must be between 1 and degree.", call. = FALSE)
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu > 1)
    stop("bgam: nu must be in (0, 1].", call. = FALSE)
  if (!is.null(mstop) && (mstop < 1L))
    stop("bgam: mstop must be >= 1.", call. = FALSE)
  if (!is.numeric(df_target) || length(df_target) != 1L || df_target <= 1)
    stop("bgam: df_target must be > 1 (target df per base-learner).",
         call. = FALSE)
  n_boot <- as.integer(n.boot)
  if (is.na(n_boot) || n_boot < 0L)
    stop("bgam: n.boot must be >= 0.", call. = FALSE)

  if (identical(family, "binomial")) {
    mean_y <- mean(y)
    if (mean_y == 0 || mean_y == 1)
      stop("bgam: response is a single class; logistic boost is undefined.",
           call. = FALSE)
  }

  # ---- Column names -------------------------------------------------------
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(p_raw))
  }
  pred_names <- colnames(x)
  p <- p_raw

  # ---- B-spline basis construction ----------------------------------------
  # For each predictor j, build:
  #   - knots_j : full knot vector (nknots + 2*degree values)
  #   - B_j     : n x K_j design matrix
  #   - K_j     : number of B-spline columns
  #   - range_j : c(x_min, x_max) for predict()
  B_list   <- vector("list", p)
  Bt_list  <- vector("list", p)
  DtD_list <- vector("list", p)
  knots_list <- vector("list", p)
  range_list <- vector("list", p)
  Kj_vec   <- integer(p)
  is_linear <- logical(p)   # near-constant predictors switched to linear

  for (j in seq_len(p)) {
    xj <- x[, j]

    # Near-constant predictor: switch to unpenalised linear base-learner
    is_near_const <- var(xj) < 1e-10 * max(vapply(seq_len(p), function(jj)
      var(x[, jj]), numeric(1L)))
    if (!is.null(unpenalized) && pred_names[j] %in% unpenalized)
      is_near_const <- TRUE

    if (is_near_const) {
      warning("bgam: predictor ", j, " (", pred_names[j],
              ") is near-constant; switching to unpenalized linear ",
              "base-learner.", call. = FALSE)
      is_linear[j] <- TRUE
      # Linear base-learner: single column (xj centred), no penalty
      xj_c <- xj - mean(xj)
      Bj <- matrix(xj_c, ncol = 1L)
      B_list[[j]]   <- Bj
      Bt_list[[j]]  <- t(Bj)
      DtD_list[[j]] <- matrix(0, 1L, 1L)
      knots_list[[j]] <- c(min(xj), max(xj))
      range_list[[j]] <- c(min(xj), max(xj))
      Kj_vec[j] <- 1L
      next
    }

    x_min <- min(xj)
    x_max <- max(xj)
    range_list[[j]] <- c(x_min, x_max)

    # Interior knots at equally-spaced quantiles
    probs    <- seq(0, 1, length.out = nknots + 2L)[seq(2L, nknots + 1L)]
    int_knots <- as.numeric(stats::quantile(xj, probs = probs))
    # Full knot vector: boundary repeats
    full_knots <- c(rep(x_min, degree), int_knots, rep(x_max, degree))
    knots_list[[j]] <- full_knots

    res    <- bgam_bspline_basis_cpp(xj, full_knots, degree)
    Bj     <- res$B
    B_list[[j]]  <- Bj
    Bt_list[[j]] <- t(Bj)
    Kj_vec[j]    <- res$K

    # Penalty matrix
    DtD_list[[j]] <- bgam_diff_penalty_cpp(res$K, dpen)
  }

  # ---- Lambda calibration -------------------------------------------------
  lambda_vec <- numeric(p)

  if (identical(lambda_method, "fixed")) {
    if (is.null(lambda_fixed))
      stop("bgam: lambda_fixed must be supplied when lambda_method='fixed'.",
           call. = FALSE)
    if (length(lambda_fixed) == 1L) {
      lambda_vec <- rep(lambda_fixed, p)
    } else if (length(lambda_fixed) == p) {
      lambda_vec <- as.numeric(lambda_fixed)
    } else {
      stop("bgam: lambda_fixed must be a scalar or a vector of length p.",
           call. = FALSE)
    }
  } else {
    # lambda_method == "df": calibrate from df_target
    for (j in seq_len(p)) {
      if (is_linear[j]) {
        lambda_vec[j] <- 0
        next
      }
      lambda_vec[j] <- .bgam_lambda_from_df(
        B_list[[j]], DtD_list[[j]], df_target
      )
    }
  }

  # ---- Cholesky prefactors (gaussian path; binomial refactors per step) ---
  # For weighted gaussian fits, include prior weights in A_j = B_j'WB_j + lambda P_j.
  chol_list <- vector("list", p)
  for (j in seq_len(p)) {
    res_c <- if (!is.null(weights)) {
      bgam_prefactor_cpp(B_list[[j]], DtD_list[[j]], lambda_vec[j], weights)
    } else {
      bgam_prefactor_cpp(B_list[[j]], DtD_list[[j]], lambda_vec[j])
    }
    chol_list[[j]] <- res_c$chol
  }

  # ---- F0 initialization --------------------------------------------------
  if (isTRUE(intercept)) {
    if (identical(family, "gaussian")) {
      F0_scalar <- mean(y)
    } else {
      my <- mean(y)
      my <- max(0.001, min(0.999, my))
      F0_scalar <- log(my / (1 - my))
    }
  } else {
    F0_scalar <- 0
  }
  F0_vec <- rep(F0_scalar, n)

  # ---- Autotune CV mstop --------------------------------------------------
  cv_result <- NULL
  if (isTRUE(autotune)) {
    cv_result <- .bgam_cv_mstop(
      x = x, y = y, family = family,
      nu = nu, mstop_max = mstop_max, nfold = nfold, seed.cv = seed.cv,
      B_list = B_list, Bt_list = Bt_list, chol_list = chol_list,
      DtD_list = DtD_list, lambda_vec = lambda_vec,
      F0_scalar = F0_scalar, nthreads = nthreads, weights = weights,
      degree = degree
    )
    mstop_opt <- cv_result$mstop_opt
    mstop_use <- mstop_opt
  } else {
    mstop_opt <- NULL
    mstop_use <- if (!is.null(mstop)) as.integer(mstop) else mstop_max
  }

  # ---- Full-data boosting loop --------------------------------------------
  fam_int <- if (identical(family, "binomial")) 1L else 0L
  boost_args <- list(
    y          = y,
    F0         = F0_vec,
    B_list     = B_list,
    Bt_list    = Bt_list,
    chol_list  = chol_list,
    DtD_list   = DtD_list,
    lambda_vec = lambda_vec,
    nu         = nu,
    mstop      = mstop_use,
    family     = fam_int,
    nthreads   = nthreads
  )
  if (!is.null(weights)) boost_args$w_prior <- weights
  boost_res <- do.call(bgam_boost_cpp, boost_args)

  # Convert C++ 0-based indices to R 1-based predictor indices.
  sel_path <- as.integer(boost_res$selection_path) + 1L
  eta_train <- as.numeric(boost_res$fitted)   # link scale

  if (identical(family, "gaussian")) {
    fitted_vals <- eta_train
    resids <- y - fitted_vals
    sigma2_hat <- if (length(resids) > p) {
      sum(resids^2) / max(1L, n - p)
    } else {
      var(resids)
    }
  } else {
    fitted_vals <- 1 / (1 + exp(-eta_train))
    # Deviance residuals
    mu_c <- pmax(1e-10, pmin(1 - 1e-10, fitted_vals))
    dev_unit <- ifelse(
      y > 0.5,
      -2 * log(mu_c),
      -2 * log(1 - mu_c)
    )
    resids <- sign(y - fitted_vals) * sqrt(pmax(dev_unit, 0))
    sigma2_hat <- NULL
  }

  # ---- Selection frequency ------------------------------------------------
  # sel_path is already 1-based; tabulate directly.
  sel_freq <- tabulate(sel_path, nbins = p) / mstop_use
  names(sel_freq) <- pred_names

  # ---- Accumulated beta list ----------------------------------------------
  beta_list_r <- boost_res$beta_list
  names(beta_list_r) <- pred_names

  # ---- base_learners info -------------------------------------------------
  base_learners <- vector("list", p)
  for (j in seq_len(p)) {
    chol_j <- if (!is.null(weights)) {
      bgam_prefactor_cpp(B_list[[j]], DtD_list[[j]], lambda_vec[j], weights)
    } else {
      bgam_prefactor_cpp(B_list[[j]], DtD_list[[j]], lambda_vec[j])
    }
    base_learners[[j]] <- list(
      knots   = knots_list[[j]],
      K       = Kj_vec[j],
      degree  = degree,
      dpen    = dpen,
      lambda  = lambda_vec[j],
      range   = range_list[[j]],
      chol    = chol_j$chol,
      # Training B matrix stored for se.fit on NULL newdata path.
      B_train = B_list[[j]]
    )
  }
  names(base_learners) <- pred_names

  # ---- Assemble fit object ------------------------------------------------
  out <- list(
    coefficients        = beta_list_r,
    fitted.values       = fitted_vals,
    linear.predictors   = eta_train,
    residuals           = resids,
    selection_path      = sel_path,
    selection_frequency = sel_freq,
    loss_path           = as.numeric(boost_res$loss_path),
    nu                  = nu,
    mstop               = mstop_use,
    mstop_opt           = mstop_opt,
    family              = family,
    base_learners       = base_learners,
    nknots              = nknots,
    degree              = degree,
    dpen                = dpen,
    df_target           = df_target,
    weights             = weights,
    n                   = n,
    p                   = p,
    predictor_names     = pred_names,
    intercept_value     = F0_scalar,
    sigma2              = sigma2_hat,
    cv                  = cv_result,
    boot                = NULL,
    call                = cl,
    terms               = NULL,
    xlevels             = NULL,
    na.action           = na_idx,
    x                   = x   # stored for bagged-predict NULL path (like ares)
  )

  # ---- Bagging ------------------------------------------------------------
  if (n_boot > 0L) {
    out <- .bgam_bagged_fit(
      out, x = x, y = y, family = family, weights = weights,
      n_boot = n_boot, seed = seed, nknots = nknots, degree = degree,
      dpen = dpen, df_target = df_target, lambda_method = lambda_method,
      lambda_fixed = lambda_fixed, nu = nu, mstop_opt = mstop_use,
      nthreads = nthreads
    )
  }

  class(out) <- "bgam"
  out
}

# ============================================================================
#  .bgam_lambda_from_df
#
# Calibrates lambda_j such that tr(S_j(lambda)) = df_target, where
# S_j(lambda) = B_j (B_j'B_j + lambda P_j)^{-1} B_j' is the hat matrix.
#
# For df_target >= K_j: returns 0 (no penalty; unconstrained base-learner).
# On uniroot failure: warns and returns lambda = 1.
# ============================================================================

.bgam_lambda_from_df <- function(Bj, Pj, df_target) {
  K <- ncol(Bj)
  if (df_target >= K) return(0)

  BtB <- crossprod(Bj)

  # tr(S(lambda)) = tr(B (B'B + lambda P)^{-1} B')
  #               = tr((B'B + lambda P)^{-1} B'B)
  #               = sum(diag(A_inv %*% BtB))
  df_fn <- function(lam) {
    A_inv <- tryCatch(
      solve(BtB + lam * Pj),
      error = function(e) NULL
    )
    if (is.null(A_inv)) return(-Inf)
    sum(diag(A_inv %*% BtB)) - df_target
  }

  # Check if df_target is achievable (avoid infinite lambda)
  # The minimum df is the penalty null space dimension (= dpen).
  # If df_fn at a very large lambda is still >= 0, set lambda to very large.
  f_lo <- df_fn(1e-6)
  f_hi <- df_fn(1e8)
  if (f_lo <= 0) return(1e-6)  # even tiny lambda gives <= df_target
  if (f_hi >= 0) return(1e8)   # even huge lambda gives >= df_target

  # f_lo > 0 and f_hi < 0: root exists in (1e-6, 1e8)
  lam <- tryCatch(
    stats::uniroot(df_fn, interval = c(1e-6, 1e8), tol = 1e-8)$root,
    error = function(e) NULL
  )
  if (is.null(lam)) {
    warning("bgam: lambda calibration (uniroot) failed for a predictor; ",
            "using lambda = 1.", call. = FALSE)
    return(1)
  }
  lam
}

# ============================================================================
#  .bgam_cv_mstop
#
# k-fold cross-validation for mstop selection.
# Grid: seq(10, mstop_max, by = 10).
# Gaussian folds: random; Binomial: stratified by y.
# Successive-halving early stop: if best CV loss at grid point t exceeds
# best_so_far by > 2% of the range of the cv_loss vector, stop early.
# ============================================================================

.bgam_cv_mstop <- function(x, y, family, nu, mstop_max, nfold, seed.cv,
                            B_list, Bt_list, chol_list, DtD_list,
                            lambda_vec, F0_scalar, nthreads, weights,
                            degree) {
  n <- nrow(x)
  p <- length(B_list)
  fam_int <- if (identical(family, "binomial")) 1L else 0L

  # Reduce mstop_max for small samples
  if (n < 100L) mstop_max <- min(mstop_max, max(50L, 2L * n))
  grid <- seq(10L, mstop_max, by = 10L)
  ng   <- length(grid)

  if (!is.null(seed.cv)) set.seed(as.integer(seed.cv))

  # Build fold assignments
  if (identical(family, "binomial")) {
    # Stratified: assign folds within each class
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)
    fold_id <- integer(n)
    # Shuffle within each stratum then assign round-robin
    ord0 <- sample(idx0)
    ord1 <- sample(idx1)
    fold_id[ord0] <- rep_len(seq_len(nfold), length(ord0))
    fold_id[ord1] <- rep_len(seq_len(nfold), length(ord1))
  } else {
    ord <- sample.int(n)
    fold_id <- integer(n)
    fold_id[ord] <- rep_len(seq_len(nfold), n)
  }

  cv_loss_mat <- matrix(NA_real_, nrow = nfold, ncol = ng)

  for (k in seq_len(nfold)) {
    te <- which(fold_id == k)
    tr <- which(fold_id != k)
    if (length(tr) < 5L || length(te) < 1L) next

    xtr <- x[tr, , drop = FALSE]
    yte <- y[te]
    wtr <- if (is.null(weights)) NULL else weights[tr]

    # Reconstruct B matrices for training fold
    B_tr  <- lapply(B_list, function(Bj) Bj[tr, , drop = FALSE])
    Bt_tr <- lapply(B_tr, t)
    DtD_l <- DtD_list

    # Recompute Cholesky for training fold
    chol_tr <- vector("list", p)
    for (j in seq_len(p)) {
      rc <- bgam_prefactor_cpp(B_tr[[j]], DtD_l[[j]], lambda_vec[j])
      chol_tr[[j]] <- rc$chol
    }

    # B matrices for test fold
    B_te <- lapply(B_list, function(Bj) Bj[te, , drop = FALSE])

    # Fit with mstop_max on training fold
    F0_tr <- rep(F0_scalar, length(tr))
    ytr <- y[tr]

    boost_k <- bgam_boost_cpp(
      y          = ytr,
      F0         = F0_tr,
      B_list     = B_tr,
      Bt_list    = Bt_tr,
      chol_list  = chol_tr,
      DtD_list   = DtD_l,
      lambda_vec = lambda_vec,
      nu         = nu,
      mstop      = mstop_max,
      family     = fam_int,
      nthreads   = nthreads
    )

    sel_k    <- as.integer(boost_k$selection_path)
    beta_k   <- boost_k$beta_list
    # beta_k holds nu-scaled accumulated betas after all mstop_max steps.
    # We reconstruct the prediction at each grid point by replaying the
    # selection path from the training fold, building up F step by step.
    # Efficient approach: replay from sel_k and accumulated increments.

    # We need nu * beta at each step separately. We rebuild from scratch.
    # Beta per step = nu * raw_step_beta. We re-run a mini-loop in R to
    # compute per-step increments efficiently.
    # The C++ engine does not return per-step betas, but we can compute
    # OOF predictions at each grid point by re-running the boost on
    # training and predicting at each grid point.

    # Since bgam_boost_cpp returns the accumulated beta after all steps,
    # we need a different approach for OOF evaluation. We run the full
    # training boost and predict the test fold using incremental reconstruction.

    # Strategy: re-run a lightweight R-level loop replaying the selection
    # path to get OOF predictions at each grid point.
    acc_beta_k <- lapply(seq_len(p), function(j) numeric(ncol(B_tr[[j]])))

    # OOF link-scale predictor for test rows: start at F0_scalar
    F_te_link <- rep(F0_scalar, length(te))

    step_idx <- 1L
    for (gi in seq_len(ng)) {
      g_target <- grid[gi]
      # Advance from current step_idx to g_target
      while (step_idx <= g_target && step_idx <= mstop_max) {
        jstar <- sel_k[step_idx] + 1L   # 1-based
        # Retrieve the step-specific beta: we need nu * beta at this step.
        # Since the C++ engine accumulates, we cannot directly retrieve
        # per-step betas from the returned beta_list. However, we can
        # use the prediction path: the contribution of step m to F_te is
        # nu * B_te[jstar] * delta_beta_m. We don't have delta_beta_m.
        #
        # Alternative efficient approach: re-run a mini C++ boost for
        # 1 step per predictor to get per-step betas. But that is too
        # expensive in the CV loop.
        #
        # Practical solution: run the full boost with mstop_max, then for
        # each grid point t, run with mstop = t and predict. Since we
        # need multiple grid points, we run the loop cumulatively in R.
        #
        # We call bgam_boost_cpp with mstop = 1 increments -- but that
        # defeats the purpose. Instead use the training beta for the full
        # run and scale by the fraction of steps.
        #
        # Correct approach (spec Section 2.6): for each grid point t,
        # run the full mstop_max training boost and then predict on the
        # test fold using only the first t steps. This requires knowing the
        # per-step betas.
        #
        # We implement this by running bgam_boost_cpp with mstop = g_target
        # (the current grid point) for each grid point. That is ng runs per
        # fold, which is at most 30 runs. This is acceptable for CV.
        step_idx <- step_idx + 1L
      }

      # Run the training fold with mstop = g_target
      boost_g <- bgam_boost_cpp(
        y          = ytr,
        F0         = F0_tr,
        B_list     = B_tr,
        Bt_list    = Bt_tr,
        chol_list  = chol_tr,
        DtD_list   = DtD_l,
        lambda_vec = lambda_vec,
        nu         = nu,
        mstop      = grid[gi],
        family     = fam_int,
        nthreads   = nthreads
      )

      # Predict on test fold
      beta_g  <- boost_g$beta_list
      F0_scl  <- F0_scalar
      F_te_g  <- bgam_predict_cpp(B_te, beta_g, F0_scl)

      # OOF loss
      if (identical(family, "gaussian")) {
        cv_loss_mat[k, gi] <- mean((yte - F_te_g)^2)
      } else {
        mu_te <- 1 / (1 + exp(-F_te_g))
        mu_te <- pmax(1e-10, pmin(1 - 1e-10, mu_te))
        cv_loss_mat[k, gi] <- -2 * mean(yte * log(mu_te) +
                                        (1 - yte) * log(1 - mu_te))
      }
    }
  }

  # Mean CV loss per grid point
  cv_loss <- colMeans(cv_loss_mat, na.rm = TRUE)

  # Successive halving early stop (re-evaluate: cap at last valid point)
  rng <- max(cv_loss, na.rm = TRUE) - min(cv_loss, na.rm = TRUE)
  best_so_far <- Inf
  last_valid <- ng
  for (gi in seq_len(ng)) {
    if (is.finite(cv_loss[gi])) {
      if (cv_loss[gi] < best_so_far) {
        best_so_far <- cv_loss[gi]
      } else if (rng > 0 && (cv_loss[gi] - best_so_far) > 0.02 * rng) {
        last_valid <- gi
        break
      }
    }
  }
  cv_loss_trimmed <- cv_loss[seq_len(last_valid)]
  grid_trimmed    <- grid[seq_len(last_valid)]
  mstop_opt       <- grid_trimmed[which.min(cv_loss_trimmed)]

  list(
    grid       = grid_trimmed,
    cv_loss    = cv_loss_trimmed,
    mstop_opt  = mstop_opt
  )
}

# ============================================================================
#  .bgam_bagged_fit
#
# Runs n.boot Poisson(1) bootstrap-weight refits with frozen hyperparameters.
# Returns the input fit object with $boot appended.
# ============================================================================

.bgam_bagged_fit <- function(fit, x, y, family, weights, n_boot, seed,
                              nknots, degree, dpen, df_target, lambda_method,
                              lambda_fixed, nu, mstop_opt, nthreads) {
  n <- nrow(x)
  w_base <- if (is.null(weights)) rep(1.0, n) else as.numeric(weights)

  has_seed <- !is.null(seed)
  if (has_seed) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(as.integer(seed))
  }

  # Use bootstrap resampling with replacement (row indices), then derive
  # weights as observation counts scaled by w_base. This matches the ares/ols
  # convention and avoids Poisson(1) weight draws that produce zeros (which
  # would be rejected by the strictly-positive weights validator).
  boot_fits <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    idx_b <- sample.int(n, n, replace = TRUE)
    # Count how many times each observation is sampled; multiply by base weights.
    counts_b <- tabulate(idx_b, nbins = n)
    w_b <- as.numeric(counts_b) * w_base
    # All zero-count obs are fine (weight = 0 means not sampled); but we
    # cannot pass zero weights to bgam.default (validator rejects w <= 0).
    # Instead, subsample to the unique rows that were drawn, with their counts.
    active <- which(counts_b > 0L)
    if (length(active) < 5L) next
    if (identical(family, "binomial") && length(unique(y[active])) < 2L) next

    f_b <- tryCatch(
      bgam.default(
        x = x[active, , drop = FALSE], y = y[active], family = family,
        nu = nu, mstop = mstop_opt, autotune = FALSE,
        nknots = nknots, degree = degree, dpen = dpen,
        df_target = df_target, lambda_method = lambda_method,
        lambda_fixed = lambda_fixed,
        weights = w_b[active], n.boot = 0L,
        nthreads = nthreads
      ),
      error = function(e) NULL
    )
    if (!is.null(f_b)) boot_fits[[b]] <- f_b
  }

  # Keep only the non-NULL fits (failed replicates are dropped).
  surviving <- Filter(Negate(is.null), boot_fits)
  fit$boot <- list(n.boot = n_boot, fits = surviving)

  # Update fitted.values and linear.predictors to reflect the bag mean.
  # Each boot fit is evaluated at the training x to get n-length predictions.
  if (length(surviving) > 0L) {
    central_preds <- as.numeric(predict.bgam(fit, newdata = x, type = "link"))
    boot_preds <- lapply(surviving, function(f) {
      tryCatch(
        as.numeric(predict.bgam(f, newdata = x, type = "link")),
        error = function(e) rep(NA_real_, nrow(x))
      )
    })
    all_link_preds <- c(list(central_preds), boot_preds)
    link_mat <- do.call(cbind, all_link_preds)
    mean_link <- rowMeans(link_mat, na.rm = TRUE)
    fit$linear.predictors <- mean_link
    if (identical(fit$family, "binomial")) {
      fit$fitted.values <- 1 / (1 + exp(-mean_link))
    } else {
      fit$fitted.values <- mean_link
    }
  }

  fit
}
