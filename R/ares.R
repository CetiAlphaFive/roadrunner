# ares -- main fitting function
#
# File layout:
#   1. ares()             -- S3 generic
#   2. ares.formula()     -- formula -> matrix dispatch
#   3. ares.default()     -- main entry; argument validation, dispatch into
#                            the C++ engine, CV / autotune / bagging
#                            orchestration, post-processing.
#   4. .term_labels()     -- earth-compatible "h(x-cut)" labels for dirs/cuts.
#   5. .ares_cv_fit()     -- CV pruning loop (pmethod = "cv").
#   6. .ares_autotune()   -- inner-CV grid search over (degree, penalty, nk,
#                            fast.k); supports warm-start subsample probe and
#                            a shared-forward fast path across grid cells.
#
# Conventions:
#   - Functions prefixed with `.` are internal helpers (not exported).
#   - The C++ engine is reached via three Rcpp exports: mars_fit_cpp,
#     mars_basis_cpp, mars_backward_only_cpp (see src/ares.cpp).
#   - Determinism contract: at a fixed `seed.cv`, fits are byte-identical
#     across `nthreads`. Any change that breaks this invariant is a bug.

#' Fast Multivariate Adaptive Regression Splines
#'
#' Fits a Multivariate Adaptive Regression Splines (MARS) model
#' (Friedman 1991). A forward pass adds pairs of hinge basis functions
#' of the form `max(0, +/-(x - knot))`; a backward pass prunes terms by
#' minimising the GCV criterion (or a K-fold CV criterion if requested).
#'
#' Two interfaces are provided. The formula method takes a model formula
#' and a data frame. The default method takes a numeric predictor matrix
#' `x` and a numeric response vector `y`.
#'
#' Fits are deterministic across thread counts: at a fixed `seed.cv`,
#' results are bit-for-bit identical regardless of `nthreads`.
#'
#' @param x A numeric matrix or data frame of predictors, or a model
#'   formula when calling the formula method.
#' @param y A numeric response vector with length `nrow(x)`. For
#'   `family = "binomial"`, a 0/1 numeric, logical, or 2-level factor.
#' @param data A data frame. Used only by the formula method.
#' @param degree Maximum interaction degree. Default `1` (additive).
#'   Use `2` or `3` for two- or three-way interactions.
#' @param nk Maximum number of basis terms in the forward pass. Default
#'   scales with `ncol(x)`: `min(200, max(20, 2 * ncol(x))) + 1`.
#' @param penalty GCV penalty per knot. Default `2` for `degree = 1`,
#'   `3` otherwise. Larger values produce sparser fits.
#' @param thresh Forward-pass early-stop threshold on relative RSS
#'   improvement. Default `0.001`.
#' @param minspan Minimum gap between knots along a variable. `0`
#'   (default) picks an automatic value from `n` and `p`.
#' @param endspan Distance from the data ends within which knots are
#'   forbidden. `0` (default) picks an automatic value.
#' @param adjust.endspan Multiplier applied to `endspan` when the
#'   candidate hinge would deepen an existing interaction. Default `1`
#'   for `gaussian`, `poisson`, `gamma`; `2` for `binomial`. Pass an
#'   explicit value to override.
#' @param auto.linpreds If `TRUE`, hinge pairs whose best knot sits at a
#'   variable's range boundary are replaced by a linear term. Default
#'   `FALSE` for `gaussian`, `poisson`, `gamma`; `TRUE` for `binomial`.
#'   Pass an explicit value to override.
#' @param fast.k Size of the Fast-MARS candidate cache (Friedman 1993).
#'   Larger values rescore more candidates each step (slower, slightly
#'   more accurate); `0` rescores every candidate every step. Default
#'   `10`.
#' @param fast.beta Age penalty for stale entries in the Fast-MARS
#'   cache. Default `1.0`. Only relevant when `fast.k > 0`.
#' @param nprune Maximum number of terms after backward pruning.
#'   Default `nk`.
#' @param pmethod Pruning method.
#'   - `"backward"` (default): minimise GCV over backward-elimination subsets.
#'   - `"none"`: keep all forward-pass terms.
#'   - `"cv"`: K-fold cross-validated subset selection. Requires `nfold > 0`
#'     (sets `nfold = 5` if not specified).
#' @param nfold Number of CV folds for `pmethod = "cv"`. Default `0`
#'   (no CV; GCV-based pruning). When `nfold > 0`, `pmethod` is
#'   promoted to `"cv"`.
#' @param ncross Number of CV repetitions (each builds a fresh fold
#'   partition; per-size mean MSE is averaged). Default `1`.
#' @param stratify If `TRUE` (default), CV folds are quantile-stratified
#'   on `y`. Useful for small `n` or skewed responses.
#' @param seed.cv Optional integer seed for the CV fold partition. Pass
#'   an integer for reproducible CV; `NULL` (default) uses the current
#'   RNG state.
#' @param cv.1se If `TRUE`, applies the one-standard-error rule under
#'   `pmethod = "cv"`: among sizes within one SE of the minimum mean
#'   CV-MSE, pick the smallest. Default `FALSE` (pick the argmin).
#' @param autotune If `TRUE`, runs an inner cross-validated grid search
#'   over `(degree, penalty, nk, fast.k)` and refits the winner on the
#'   full data. The chosen settings and full grid scores are returned in
#'   `$autotune`. Default `FALSE`.
#' @param autotune.speed Speed/quality trade-off for `autotune`. Only
#'   used when `autotune = TRUE`.
#'   - `"balanced"` (default): explores a moderate `fast.k` range and
#'     picks the smallest whose CV-MSE is within 1% of the best.
#'   - `"quality"`: disables the Fast-MARS cache (most thorough,
#'     slowest).
#'   - `"fast"`: forces an aggressive cache (cheapest, slight accuracy
#'     trade-off).
#' @param autotune.warmstart If `TRUE` (default), `autotune` first tunes
#'   on a small subsample. If one cell wins decisively, the full grid is
#'   skipped and the winner is refit on all rows. Only applies when
#'   `autotune = TRUE` and `n >= 200`.
#' @param n.boot Number of bootstrap replicate fits for bagging.
#'   Default `0` (no bagging). When `> 0`, `predict()` averages over the
#'   replicates plus the central fit and (with `se.fit = TRUE`) returns
#'   a per-prediction bag standard deviation. Composes with `autotune`:
#'   each replicate reuses the central fit's chosen hyperparameters.
#' @param na.action Strategy for missing values in `x`. `"impute"`
#'   (default) replaces each column's `NA`s with that column's training
#'   median and stores the medians for `predict()` to reuse. `"omit"`
#'   drops rows with any `NA`. Either action warns. Missing values in
#'   `y`, or `NaN`/`Inf` anywhere in `x`, are always rejected.
#' @param family Response family for the final coefficient fit.
#'   - `"gaussian"` (default): identity link, OLS.
#'   - `"binomial"`: logit link. `y` may be 0/1 numeric, logical, or a
#'     2-level factor.
#'   - `"poisson"`: log link. `y` must be non-negative integer-valued.
#'   - `"gamma"`: log link. `y` must be strictly positive.
#'   For non-gaussian families, term selection runs on the numeric `y`
#'   (fast); the selected basis is then refit on the response scale via
#'   `stats::glm.fit()` with the requested family.
#' @param weights Optional non-negative numeric vector of length
#'   `nrow(x)` for weighted least-squares fitting. Default `NULL`
#'   (unweighted). Term selection, GCV / CV pruning, and autotune all
#'   respect the weights. Negative, `NA`, or `NaN` weights are rejected.
#' @param varmod Residual variance model used by `predict()` for
#'   prediction intervals (gaussian only).
#'   - `"none"` (default): no variance model is stored; `interval = "pint"`
#'     is unavailable.
#'   - `"const"`: stores a single residual SD; intervals are
#'     `yhat +/- qt(level, df) * sigma`.
#'   - `"lm"`: fits a small linear model of `|resid|` on `yhat` to allow
#'     simple **yhat-dependent** heteroscedasticity. Captures residual scale
#'     that changes linearly with the fitted mean. Does NOT capture residual
#'     scale that depends on a predictor whose contribution to `yhat` is
#'     small -- in that case the |resid| ~ yhat slope is close to zero and
#'     `"lm"` collapses to roughly the same as `"const"`. If you suspect
#'     x-driven heteroscedasticity (e.g. variance depends on a covariate
#'     orthogonal to the mean structure), `varmod = "lm"` will not help
#'     and coverage will degrade in high-variance regions.
#'   Ignored for non-gaussian families.
#' @param trace Trace level. `0` (default) is silent; `1` reports
#'   forward-pass progress.
#' @param nthreads Number of threads. `0` (default) uses
#'   `RcppParallel::defaultNumThreads()`. Examples and the vignette cap
#'   this at `2` for CRAN compliance.
#' @param ... Currently ignored.
#' @return An object of class `"ares"`: a list containing the fitted
#'   coefficients, the basis matrix `bx`, the term directions `dirs`
#'   and knots `cuts`, the indices of `selected.terms`, training `rss`
#'   and `gcv`, `fitted.values` and `residuals`, predictor names
#'   `namesx`, the call, and echoed control parameters. Non-gaussian
#'   fits additionally include `$family`, `$glm` (with `deviance`,
#'   `null.deviance`, `df.null`, `df.residual`, `aic`, `converged`,
#'   `iter`), and `$linear.predictor`; `$fitted.values` is on the
#'   response scale (probabilities for binomial, positive means for
#'   poisson and gamma). Autotune fits carry `$autotune` and bagged
#'   fits carry `$boot`.
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics* 19(1):1-67.
#'
#' Friedman, J. H. (1993). *Fast MARS*. Stanford University Department
#' of Statistics Technical Report 110.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' print(fit)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
# ============================================================================
#  S3 generic + formula method
# ============================================================================

#' @export
ares <- function(x, ...) UseMethod("ares")

#' @rdname ares
#' @export
ares.formula <- function(x, data = NULL, ..., y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()
  # Use `na.action = na.pass` so model.frame doesn't blow up on NAs; the
  # default method then applies its own (median-impute or omit) rule per
  # its `na.action` argument. y-side NAs still trigger a hard stop in
  # ares.default(), which is correct -- model.response should fail loud
  # when the target is missing.
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  yv <- stats::model.response(mf)
  if (is.null(yv)) stop("ares: response variable is missing from formula/data.")
  mm <- stats::model.matrix(formula, mf)
  # drop intercept column
  has_int <- "(Intercept)" %in% colnames(mm)
  if (has_int) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  # When the caller passed family = "binomial" we must keep the response in
  # its original type (factor / 0-1 / logical) so ares.default can validate
  # and coerce it consistently. For gaussian (default), preserve the old
  # behaviour of coercing to numeric.
  dots <- list(...)
  fam_in <- dots$family
  if (!is.null(fam_in) && identical(fam_in, "binomial")) {
    out <- ares.default(x = mm, y = yv, ...)
  } else {
    out <- ares.default(x = mm, y = as.numeric(yv), ...)
  }
  out$call <- cl
  out$terms <- stats::terms(formula, data = data)
  out
}

# ============================================================================
#  Main default method
# ============================================================================
#
# Orchestrates: input validation, default selection, optional autotune,
# optional CV pruning, optional bagging, and the C++ engine call(s).
#
# Returns an object of class "ares" with the components documented at the
# top of the file.

#' @rdname ares
#' @export
ares.default <- function(x, y, degree = 1L, nk = NULL, penalty = NULL,
                         thresh = 0.001, minspan = 0L, endspan = 0L,
                         adjust.endspan = 1L, auto.linpreds = FALSE,
                         fast.k = 10L, fast.beta = 1.0,
                         nprune = NULL,
                         pmethod = c("backward", "none", "cv"),
                         nfold = 0L, ncross = 1L, stratify = TRUE,
                         seed.cv = NULL, cv.1se = FALSE,
                         autotune = FALSE,
                         autotune.speed = c("balanced", "quality", "fast"),
                         autotune.warmstart = TRUE,
                         n.boot = 0L,
                         na.action = c("impute", "omit"),
                         family = c("gaussian", "binomial",
                                    "poisson", "gamma"),
                         weights = NULL,
                         varmod = c("none", "const", "lm"),
                         trace = 0L, nthreads = 0L, ...) {
  cl <- match.call()
  pmethod <- match.arg(pmethod)
  autotune_speed <- match.arg(autotune.speed)
  na.action <- match.arg(na.action)
  family <- match.arg(family)
  varmod <- match.arg(varmod)
  # ---- Response coercion for `family` --------------------------------------
  # For `family = "binomial"`, we need y to be a 0/1 numeric on entry to the
  # downstream gaussian engine. Accept:
  #   - integer / numeric vectors whose unique non-NA values are subset of
  #     {0, 1}.
  #   - 2-level factors. Reference level (level 1) maps to 0; the other to
  #     1. The level pair is stashed for predict.ares() so that predicted
  #     factor responses round-trip to the original labels.
  # Anything else is rejected with a clear stop().
  y_levels <- NULL
  if (family == "binomial") {
    if (is.logical(y)) {
      y <- as.integer(y)
    } else if (is.factor(y)) {
      lv <- levels(y)
      if (length(lv) != 2L)
        stop("ares: family = 'binomial' requires y to be a 2-level factor",
             " (got ", length(lv), " levels: ",
             paste(lv, collapse = ", "), ").")
      y_levels <- lv
      y <- as.integer(y) - 1L         # 1st level -> 0, 2nd level -> 1
    } else if (is.character(y)) {
      yf <- factor(y)
      lv <- levels(yf)
      if (length(lv) != 2L)
        stop("ares: family = 'binomial' requires y to be a 2-level",
             " character vector (got ", length(lv), " levels).")
      y_levels <- lv
      y <- as.integer(yf) - 1L
    } else if (is.numeric(y) || is.integer(y)) {
      # 0/1 integers / doubles only. NaN/Inf are caught later by the
      # is.finite() check; this block only rejects out-of-{0,1} values.
      yv <- y[is.finite(y)]
      uv <- unique(yv)
      if (!all(uv %in% c(0, 1)))
        stop("ares: family = 'binomial' requires y in {0, 1}; got values ",
             paste(utils::head(sort(uv), 6), collapse = ", "), ".")
      y <- as.integer(y)
    } else {
      stop("ares: family = 'binomial' requires y to be 0/1 numeric,",
           " logical, character, or 2-level factor (got class: ",
           paste(class(y), collapse = "/"), ").")
    }
  } else if (family == "poisson") {
    # poisson: y must be nonneg, finite, integer-valued (NaN/Inf caught later
    # by the is.finite() guard). We coerce to numeric (glm.fit accepts that)
    # but enforce "integer-valued" by comparing to round(y) within a tight
    # tolerance to allow y stored as double 0.0/1.0/2.0/... .
    if (!is.numeric(y) && !is.integer(y))
      stop("ares: family = 'poisson' requires numeric y.")
    y <- as.numeric(y)
    yv <- y[is.finite(y)]
    if (any(yv < 0))
      stop("ares: family = 'poisson' requires y >= 0; got negative values.")
    if (any(abs(yv - round(yv)) > 1e-8))
      stop("ares: family = 'poisson' requires integer-valued y;",
           " got non-integer values (e.g. ",
           utils::head(yv[abs(yv - round(yv)) > 1e-8], 3)[1], ").")
  } else if (family == "gamma") {
    # gamma (log link): y must be strictly positive and finite.
    if (!is.numeric(y) && !is.integer(y))
      stop("ares: family = 'gamma' requires numeric y.")
    y <- as.numeric(y)
    yv <- y[is.finite(y)]
    if (any(yv <= 0))
      stop("ares: family = 'gamma' requires y > 0; got non-positive values.")
  }

  # ---- Family-conditional default flip (binomial only) --------------------
  # For family == "binomial" the empirical sweet spot is the earth-like
  # configuration auto.linpreds = TRUE, adjust.endspan = 2. The A/B test in
  # inst/sims/v0.25_ab_binomial.R showed this ties or beats the conservative
  # defaults on all mlbench binary cells (4 DGPs x 2 sample sizes). Gaussian
  # remains conservative — the earth-like flip regresses regression cells
  # in the v0.20-mlbench grid. We only flip arguments the user did not
  # explicitly pass; users who pin either argument retain their choice.
  if (family == "binomial") {
    arg_names <- names(cl)
    if (!("auto.linpreds" %in% arg_names)) {
      auto.linpreds <- TRUE
    }
    if (!("adjust.endspan" %in% arg_names)) {
      adjust.endspan <- 2L
    }
  }

  # ---- Observation weights (initial validation) ----------------------------
  # Validate up front (length, finiteness, non-negativity) so subsequent
  # row-subsetting branches (na.action = "omit") can subset weights too.
  # Normalisation to mean(w) = 1 happens later, after any row drops.
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("ares: weights must be numeric.")
    if (length(weights) != length(y))
      stop("ares: length(weights) (", length(weights),
           ") must equal length(y) (", length(y), ").")
    if (any(!is.finite(weights)))
      stop("ares: weights contain NA / NaN / Inf; all weights must be finite.")
    # BUG-005 (v0.0.0.9029): require strictly positive weights. Zero
    # weights used to be permitted, but the C++ engines GCV denominator
    # (and varmod df, and auto_minspan / auto_endspan) all use the raw
    # row count n -- under sum(w)=n normalisation, n_effective should be
    # sum(w_i > 0). With zero-weight rows present, n > n_effective and
    # the GCV penalty per knot is implicitly under-estimated, causing
    # over-fitting (probe 12 in the audit: 50 zero-weight rows kept 10
    # terms vs 7 in the dropped-rows fit). Rejecting up front is
    # lossless: the user can drop those rows from x/y themselves with
    # identical semantics, and the determinism invariant is preserved
    # because no C++ engine code changes.
    if (any(weights <= 0))
      stop("ares: weights must be strictly positive. To omit a row, drop",
           " it from x/y; weights = 0 is not a substitute (it would bias",
           " GCV and the variance-model df under the engines",
           " sum(w) = n normalisation).")
  }
  # nfold > 0 implies CV pruning. Promote pmethod to "cv" so downstream
  # branches share the same code path; a user explicitly asking for
  # pmethod="cv" with nfold==0 gets a default of 5 folds.
  nfold <- as.integer(nfold)
  ncross <- as.integer(ncross)
  if (is.na(nfold)  || nfold  < 0L) nfold  <- 0L
  if (is.na(ncross) || ncross < 1L) ncross <- 1L
  if (pmethod == "cv" && nfold == 0L) nfold <- 5L
  if (pmethod != "cv" && nfold >  0L) pmethod <- "cv"
  # Internal int code for the C++ engine. The CV path uses the GCV-based
  # backward pass underneath (pmethod_int = 0) and overrides the size pick
  # from the R-side via force_size; pmethod="none" keeps the old behaviour.
  pmethod_int <- if (pmethod == "none") 1L else 0L

  # ---- Input validation ----
  # When `x` is a data.frame with factor or character columns, expand it to
  # a numeric model matrix via `model.matrix(~ . , x)` (treatment contrasts,
  # k-1 dummies per factor; intercept dropped). This matches the formula
  # method's path so the two interfaces behave the same on mixed-type data.
  # Levels are stashed on `$factor_info` so predict.ares() can replay the
  # same expansion on newdata.
  factor_info <- NULL
  pre_na_medians <- NULL    # numeric-column medians captured BEFORE factor
                            # expansion (so we can replay them at predict
                            # time, before model.matrix runs there too).
  if (is.data.frame(x)) {
    is_cat <- vapply(x, \(z) is.factor(z) || is.character(z), logical(1L))
    is_num <- vapply(x, is.numeric, logical(1L))

    # ---- Pre-expansion NA handling on numeric columns -------------------
    # If we call model.matrix(~ ., x) with NAs in any numeric column, the
    # default na.action is na.omit which silently drops rows -- and the
    # `na.action` arg we honour at the matrix level below would never see
    # those rows. Handle numeric-column NAs HERE so the expansion is
    # clean. Factor / character columns with NAs are left to model.matrix
    # (which already errors loudly; rare enough that the strict path is
    # the right default).
    if (any(is_num)) {
      num_names <- names(x)[is_num]
      df_na_mask <- vapply(num_names, \(jn) anyNA(x[[jn]]), logical(1L))
      if (any(df_na_mask)) {
        rows_aff <- sum(rowSums(is.na(as.matrix(x[, num_names, drop = FALSE]))) > 0L)
        cols_aff <- sum(df_na_mask)
        if (na.action == "impute") {
          pre_na_medians <- vapply(num_names,
                                   \(jn) stats::median(x[[jn]], na.rm = TRUE),
                                   numeric(1L))
          names(pre_na_medians) <- num_names
          if (any(!is.finite(pre_na_medians)))
            stop("ares: na.action='impute' but at least one numeric column",
                 " is entirely NA; cannot compute median. Drop the column",
                 " or pass na.action='omit'.")
          for (jn in num_names) {
            na_j <- is.na(x[[jn]])
            if (any(na_j)) x[[jn]][na_j] <- pre_na_medians[[jn]]
          }
          warning("Missing values in X: median-imputed ", sum(df_na_mask),
                  " numeric column(s) across ", rows_aff, " row(s).",
                  " Column medians are stored on the fit and reapplied at",
                  " predict() time. Pass `na.action='omit'` to drop",
                  " incomplete rows instead.", call. = FALSE)
        } else { # "omit"
          keep <- stats::complete.cases(x[, num_names, drop = FALSE])
          dropped <- sum(!keep)
          warning("Missing values in X: dropped ", dropped,
                  " incomplete row(s) (across ", cols_aff,
                  " numeric column(s) with NA). Pass `na.action='impute'`",
                  " to median-impute instead.", call. = FALSE)
          x <- x[keep, , drop = FALSE]
          y <- y[keep]
          if (!is.null(weights)) weights <- weights[keep]
          if (length(y) < 3L)
            stop("ares: na.action='omit' left fewer than 3 rows; aborting.")
        }
      }
    }

    if (any(is_cat)) {
      # Coerce character columns to factor with their observed levels (the
      # natural lexical order; users wanting a specific reference level
      # should pre-cast to factor with the desired `levels`).
      for (j in which(is_cat)) {
        if (is.character(x[[j]])) x[[j]] <- factor(x[[j]])
      }
      xlevels <- lapply(x[is_cat], levels)
      xmm <- stats::model.matrix(~ ., data = x)
      if ("(Intercept)" %in% colnames(xmm))
        xmm <- xmm[, colnames(xmm) != "(Intercept)", drop = FALSE]
      factor_info <- list(
        xlevels  = xlevels,                 # named list of levels per factor col
        orig_names = names(x),              # original df column names
        is_cat   = is_cat,                  # logical(ncol(x_df))
        is_num   = is_num,                  # logical(ncol(x_df))
        expanded_names = colnames(xmm),     # expanded dummy column names
        # Save numeric medians keyed by ORIGINAL column name so predict.ares
        # can fill newdata's NAs before model.matrix runs there. (The
        # post-expansion median map in $na.medians is keyed by expanded
        # column name and applies only when the matrix path is used.)
        num_medians = pre_na_medians
      )
      x <- xmm
    } else {
      nm <- names(x)
      x <- as.matrix(x)
      storage.mode(x) <- "double"
      colnames(x) <- nm
    }
  }
  if (!is.matrix(x) || !is.numeric(x))
    stop("ares: x must be a numeric matrix or data frame.")
  if (!is.numeric(y))
    stop("ares: y must be numeric.")
  # Missing values in y are a hard error (no sensible imputation when y is
  # the regression target). NaN / +-Inf are also rejected — they would
  # propagate through OLS and break the fit.
  if (any(!is.finite(y)))
    stop("ares: y contains NA, NaN, or non-finite values; aborting fit.")
  if (length(y) != nrow(x))
    stop("ares: length(y) (", length(y), ") must equal nrow(x) (", nrow(x), ").")

  # ---- Missing-X handling --------------------------------------------------
  # NaN / +-Inf in x are still rejected outright (they break the basis
  # evaluation arithmetic regardless of any imputation rule). NA values are
  # handled by `na.action`:
  #   - "impute" (default): fill each column's NAs with that column's
  #     median (computed from the non-missing rows). Column medians are
  #     stashed on the fit as `$na.medians` so predict.ares() can apply
  #     the same imputation to future newdata.
  #   - "omit": drop every row that has any NA in x.
  # Either action emits a warning summarising what happened.
  na_medians <- NULL
  if (any(is.nan(x)) || any(is.infinite(x)))
    stop("ares: x contains NaN or +/-Inf values; aborting fit.",
         " (NA is handled via `na.action`; NaN / Inf are not.)")
  na_x <- is.na(x)
  if (any(na_x)) {
    rows_aff <- sum(rowSums(na_x) > 0L)
    cols_aff <- sum(colSums(na_x) > 0L)
    if (na.action == "impute") {
      na_medians <- vapply(seq_len(ncol(x)),
                           \(j) stats::median(x[, j], na.rm = TRUE),
                           numeric(1L))
      names(na_medians) <- colnames(x)
      if (any(!is.finite(na_medians)))
        stop("ares: na.action='impute' but at least one x column is",
             " entirely NA; cannot compute median. Drop the column or",
             " pass na.action='omit'.")
      for (j in seq_len(ncol(x))) {
        na_j <- na_x[, j]
        if (any(na_j)) x[na_j, j] <- na_medians[j]
      }
      warning("Missing values in X: median-imputed ",
              sum(na_x), " NA value(s) across ", rows_aff,
              " row(s) and ", cols_aff, " column(s). Column medians are",
              " stored on the fit and reapplied at predict() time.",
              " Pass `na.action='omit'` to drop incomplete rows instead.",
              call. = FALSE)
    } else { # "omit"
      keep <- rowSums(na_x) == 0L
      dropped <- sum(!keep)
      warning("Missing values in X: dropped ", dropped,
              " incomplete row(s) (across ", cols_aff,
              " column(s) with NA). Pass `na.action='impute'` to",
              " median-impute instead.", call. = FALSE)
      x <- x[keep, , drop = FALSE]
      y <- y[keep]
      if (!is.null(weights)) weights <- weights[keep]
      if (length(y) < 3L)
        stop("ares: na.action='omit' left fewer than 3 rows; aborting.")
    }
  }
  if (!is.numeric(degree) || length(degree) != 1L || degree < 1)
    stop("ares: degree must be a single integer >= 1.")
  degree <- as.integer(degree)
  if (ncol(x) < 1L)
    stop("ares: x must have at least one column.")

  # ---- Drop constant columns with warning ----
  col_const <- vapply(seq_len(ncol(x)), function(j) {
    z <- x[, j]
    diff(range(z)) < .Machine$double.eps * 16
  }, logical(1))
  dropped_names <- character(0)
  if (any(col_const)) {
    dropped_names <- colnames(x)[col_const]
    if (is.null(dropped_names)) dropped_names <- paste0("V", which(col_const))
    warning("ares: dropped ", sum(col_const), " constant column(s): ",
            paste(dropped_names, collapse = ", "), call. = FALSE)
    x <- x[, !col_const, drop = FALSE]
    if (ncol(x) < 1L)
      stop("ares: all predictor columns are constant after dropping.")
  }

  # ---- Observation weights (final normalisation) --------------------------
  # Normalise to mean(w) = 1 so that sum(w) == n. The C++ engine assumes this
  # convention (so its GCV denominator and intercept-init math stay valid).
  # Constant-column drops above don't change weight length; we operate on the
  # current (possibly omit-subsetted) `weights`.
  w_norm <- NULL
  if (!is.null(weights)) {
    wm <- mean(weights)
    if (wm <= 0)
      stop("ares: post-NA-omit weights have zero mean; cannot fit.")
    w_norm <- as.numeric(weights) / wm
  }

  # ---- Defaults ----
  p <- ncol(x)
  if (is.null(nk)) nk <- min(200L, max(20L, 2L * p)) + 1L
  nk <- as.integer(nk)
  if (nk < 3L) {
    warning("ares: nk coerced to 3 (minimum).", call. = FALSE)
    nk <- 3L
  }
  if (is.null(penalty)) penalty <- if (degree > 1L) 3 else 2
  penalty <- as.numeric(penalty)
  if (is.null(nprune)) nprune <- nk
  nprune <- as.integer(nprune)
  trace <- as.integer(trace)
  minspan <- as.integer(minspan)
  endspan <- as.integer(endspan)
  adjust_endspan <- as.integer(adjust.endspan)
  if (adjust_endspan < 1L) {
    warning("ares: adjust.endspan coerced to 1 (minimum).", call. = FALSE)
    adjust_endspan <- 1L
  }
  auto_linpreds <- as.integer(isTRUE(auto.linpreds))
  fast_k <- as.integer(fast.k)
  if (is.na(fast_k) || fast_k < 0L) fast_k <- 0L  # 0 = unlimited (no caching)
  fast_beta <- as.numeric(fast.beta)
  if (is.na(fast_beta) || fast_beta < 0) fast_beta <- 0
  n_boot <- as.integer(n.boot)
  if (is.na(n_boot) || n_boot < 0L) n_boot <- 0L

  nthreads <- as.integer(nthreads)
  nthreads_eff <- if (nthreads <= 0L) RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  # ---- Capture column names ----
  namesx <- colnames(x)
  if (is.null(namesx)) namesx <- paste0("V", seq_len(p))

  # ---- Coerce x to plain numeric matrix ----
  storage.mode(x) <- "double"

  # ---- Autotune dispatch (Phase 2 -- v0.15+) ----
  if (isTRUE(autotune)) {
    # Promote nfold to a sensible default if user didn't ask. v0.0.0.9023:
    # on high-p problems (p >= 15, where nk_eff >= 31), default to nfold = 3
    # instead of 5. Inner-CV variance from 3 folds vs 5 is within wash on
    # the standard mlbench DGPs (Friedman-1 / additive / interaction at
    # p = 20: MSE shift <2.5%), and per-fold cost dominates wall-clock so
    # 3 folds gives ~40% reduction on top of the v0.0.0.9021/9022 cuts.
    p_eff <- ncol(x)
    nfold_default <- if (p_eff >= 15L) 3L else 5L
    nfold_at <- if (nfold > 0L) nfold else nfold_default
    out <- .ares_autotune(
      x = x, y = as.numeric(y),
      nk = as.integer(nk),
      thresh = thresh, minspan = minspan, endspan = endspan,
      adjust_endspan = adjust_endspan, auto_linpreds = auto_linpreds,
      fast_k = fast_k, fast_beta = fast_beta,
      nprune = as.integer(nprune),
      nfold = nfold_at, ncross = ncross, stratify = isTRUE(stratify),
      seed_cv = seed.cv, cv_1se = isTRUE(cv.1se),
      trace = trace, nthreads = nthreads_eff,
      autotune_speed = autotune_speed,
      warmstart = isTRUE(autotune.warmstart),
      weights = w_norm
    )
  } else if (pmethod == "cv") {
    out <- .ares_cv_fit(x, y = as.numeric(y),
                        degree = degree, nk = as.integer(nk),
                        penalty = as.numeric(penalty), thresh = thresh,
                        minspan = minspan, endspan = endspan,
                        adjust_endspan = adjust_endspan,
                        auto_linpreds = auto_linpreds,
                        fast_k = fast_k, fast_beta = fast_beta,
                        nprune = as.integer(nprune),
                        trace = trace, nthreads = nthreads_eff,
                        nfold = nfold, ncross = ncross,
                        stratify = isTRUE(stratify),
                        seed_cv = seed.cv,
                        cv_1se = isTRUE(cv.1se),
                        weights = w_norm)
  } else {
    out <- mars_fit_cpp(x, as.numeric(y), degree, nk, penalty, thresh,
                        minspan, endspan, adjust_endspan, auto_linpreds,
                        fast_k, fast_beta,
                        nprune, pmethod_int, trace, nthreads_eff,
                        force_size = 0L, return_path = 0L,
                        weights_in = w_norm)
  }

  # ---- Post-process ----
  rownames(out$dirs) <- .term_labels(out$dirs, out$cuts, namesx)
  rownames(out$cuts) <- rownames(out$dirs)
  colnames(out$dirs) <- namesx
  colnames(out$cuts) <- namesx
  colnames(out$bx) <- rownames(out$dirs)[out$selected.terms]
  names(out$coefficients) <- colnames(out$bx)
  out$fitted.values <- drop(out$bx %*% out$coefficients)
  out$residuals <- as.numeric(y) - out$fitted.values
  out$namesx <- namesx
  out$call <- cl
  out$pmethod <- if (isTRUE(autotune)) "backward" else pmethod
  out$dropped <- dropped_names

  # ---- Bagging (Phase 2 -- v0.18+) ----
  # Refit on n_boot row-bootstrap samples using the central fit's tuned
  # hyperparameters. Stored in $boot$fits as plain mars_fit_cpp returns
  # (with dirs/cuts/selected.terms suitable for predict via mars_basis_cpp).
  if (n_boot > 0L) {
    # Use the chosen (degree, penalty, nk, fast.k) from the central fit so
    # bagging doesn't re-run autotune for every replicate.
    deg_b <- if (isTRUE(autotune)) out$autotune$degree  else degree
    pen_b <- if (isTRUE(autotune)) out$autotune$penalty else penalty
    nk_b  <- if (isTRUE(autotune)) out$autotune$nk      else nk
    fk_b  <- if (isTRUE(autotune)) out$autotune$fast_k  else fast_k
    pmethod_int_b <- if (out$pmethod == "none") 1L else 0L

    # RNG protection so bagging only consumes the user's RNG when seed.cv
    # is unset.
    has_seed <- !is.null(seed.cv)
    if (has_seed) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(),
                         inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
      }
      # Offset by a fixed prime so bagging RNG stream is distinct from CV's.
      set.seed(as.integer(seed.cv) + 1009L)
    }
    n_full <- length(y)
    boot_fits <- vector("list", n_boot)
    # BUG-001 (v0.0.0.9029): capture each bag's bootstrap indices up-front so the
    # post-hoc GLM refit (binomial / poisson / gamma) below re-uses the SAME
    # rows the basis was selected on. Previously the refit loop redrew its own
    # indices from the live RNG and got a different bootstrap sample whenever
    # `seed.cv = NULL` (the unseeded path).
    idx_list <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n_full, n_full, replace = TRUE)
      idx_list[[b]] <- idx
      x_b <- x[idx, , drop = FALSE]
      y_b <- as.numeric(y)[idx]
      # Bootstrap-resampled weights, renormalised so mean(w_b) = 1 inside the
      # replicate. NULL when no weights were supplied → bag fits behave
      # exactly like the pre-weights bagging path.
      w_b <- if (!is.null(w_norm)) {
        wb <- w_norm[idx]
        wb / mean(wb)
      } else NULL
      f_b <- mars_fit_cpp(
        x_b, y_b, as.integer(deg_b), as.integer(nk_b),
        as.numeric(pen_b),
        thresh, minspan, endspan, adjust_endspan, auto_linpreds,
        as.integer(fk_b), fast_beta,
        if (length(nprune) == 0L || is.na(nprune[1])) as.integer(nk_b) else as.integer(nprune),
        pmethod_int_b, trace, nthreads_eff, 0L, 0L,
        weights_in = w_b
      )
      # Strip heavy fields that bagging doesn't need (bx is n_full x M).
      f_b$bx <- NULL
      boot_fits[[b]] <- f_b
    }
    out$boot <- list(
      fits   = boot_fits,
      n.boot = n_boot,
      idx    = idx_list
    )
  }

  # Persist the NA-handling metadata so predict.ares() can reapply the
  # same imputation to newdata. `na.medians` is NULL when the training
  # data had no NAs or `na.action = "omit"` was used.
  out$na.action <- na.action
  out$na.medians <- na_medians
  # BUG-002 (v0.0.0.9029): keep the (post-NA-imputed, post-factor-expansion,
  # post-constant-drop) training x on the fit so `predict(fit)` -- the
  # natural form with `newdata = NULL` -- can re-route through the bag-
  # averaging path for bagged fits. Without this the NULL branch
  # short-circuited to `$fitted.values` (the central fit only), silently
  # diverging from `predict(fit, x_train)` on bagged fits. Memory cost is
  # 8 * n * p bytes; cheap relative to bx (n * M) which we already store.
  out$x <- x
  # Factor-expansion metadata: NULL when training x was already a numeric
  # matrix (or a data.frame with no factor/character columns). When set,
  # predict.ares() replays the same model.matrix expansion on newdata.
  out$factor_info <- factor_info

  # ---- family bookkeeping --------------------------------------------------
  out$family <- family
  # `w_norm` is the same vector passed to the C++ engine (mean(w_norm) = 1
  # or NULL); we pass it onward to glm.fit() so the post-hoc GLM solves the
  # weighted-MLE on the selected basis. For NULL weights, glm.fit defaults
  # to equal weights, matching the pre-weights binomial path exactly.
  if (family == "binomial") {
    # Post-hoc GLM refit on the selected basis (earth strategy). The
    # forward + backward already picked terms via OLS on the binary y; now
    # IRLS estimates logit-scale coefficients on the same design. We use
    # `glm.fit` (no formula, no model.frame) for speed.
    out <- .ares_refit_binomial(out, y_int = as.integer(y),
                                y_levels = y_levels,
                                weights = w_norm)
    # Refit each bag replicate too (each replicate already holds its own
    # selected basis on its bootstrap rows; we re-evaluate on those rows
    # and IRLS-refit). The original x/y are needed because bag fits only
    # carry dirs/cuts/selected.terms (not their bootstrap data).
    if (!is.null(out$boot) && length(out$boot$fits)) {
      out$boot$fits <- .ares_refit_boot_binomial(out$boot$fits,
                                                  x_full = x,
                                                  y_full = as.integer(y),
                                                  idx_list = out$boot$idx,
                                                  weights = w_norm)
    }
  } else if (family == "poisson" || family == "gamma") {
    # Mirror of the binomial path. Forward + backward already ran as
    # gaussian on the selected basis; now refit the chosen basis under the
    # appropriate GLM (poisson(log) or Gamma(log)) so the coefficients and
    # fitted means come out on the correct scale.
    out <- .ares_refit_glm_loglink(out, y = as.numeric(y),
                                   family_name = family,
                                   weights = w_norm)
    if (!is.null(out$boot) && length(out$boot$fits)) {
      out$boot$fits <- .ares_refit_boot_glm_loglink(out$boot$fits,
                                                    x_full = x,
                                                    y_full = as.numeric(y),
                                                    family_name = family,
                                                    idx_list = out$boot$idx,
                                                    weights = w_norm)
    }
  }

  # ---- Variance model (gaussian only) -------------------------------------
  # Used by predict.ares(interval = "pint"). varmod = "none" (default)
  # leaves the slot empty; binomial / poisson / gamma always skip this
  # because their predictive variance is mean-determined and we do not
  # implement family-specific PIs here.
  out$varmod <- NULL
  if (varmod != "none" && family == "gaussian") {
    out$varmod <- .ares_fit_varmod(out, varmod = varmod, weights = w_norm)
  }

  class(out) <- c("ares")
  out
}

# ============================================================================
#  Internal: post-hoc GLM refit (family = "binomial")
# ============================================================================
#
# Refits a binomial GLM on the *selected* basis (`out$bx`) using IRLS via
# `stats::glm.fit`. Replaces `$coefficients`, `$fitted.values`,
# `$linear.predictor`, `$residuals`, and `$glm`. The forward + backward
# components (`dirs`, `cuts`, `selected.terms`, `rss`, `gcv`) are kept
# untouched — they reflect the gaussian model-selection stage.
#
# @keywords internal
.ares_refit_binomial <- function(out, y_int, y_levels, weights = NULL) {
  bx <- out$bx
  if (is.null(bx) || nrow(bx) == 0L)
    stop("ares: family='binomial' but $bx is empty; cannot fit GLM.")
  # The intercept column is the first basis term (all ones), so glm.fit gets
  # the full design and we pass `intercept = FALSE` to avoid duplicating it.
  # Suppress glm.fit's own warnings about non-convergence / numerically 0-1
  # probabilities -- we surface those via $glm$converged so callers can act
  # on them, but in normal use (well-separated DGP at small n) they are
  # routine and shouldn't spam the console.
  g <- suppressWarnings(stats::glm.fit(
    x = bx, y = y_int, family = stats::binomial(),
    weights = if (is.null(weights)) rep(1, length(y_int)) else weights,
    intercept = FALSE,
    control = list(maxit = 50L, epsilon = 1e-8, trace = FALSE)
  ))
  cf <- g$coefficients
  # On rank-deficient basis or perfect separation, glm.fit can emit NA
  # coefficients. Replace NA with 0 so downstream prediction stays finite
  # (this is the same fallback earth uses internally; the affected terms
  # contribute nothing to the linear predictor).
  if (any(!is.finite(cf))) cf[!is.finite(cf)] <- 0
  names(cf) <- colnames(bx)
  # Clamp the linear predictor to +/- 30 so plogis stays in (1e-13, 1-1e-13).
  # This is a defensive bound; with a well-conditioned design it never
  # activates. Earth applies a similar clamp inside its predict.earth().
  lp <- drop(bx %*% cf)
  lp_clamped <- pmin(pmax(lp, -30), 30)
  out$coefficients   <- cf
  out$linear.predictor <- lp_clamped
  out$fitted.values  <- stats::plogis(lp_clamped)
  # Residuals = response - probability (working / response residuals).
  out$residuals      <- as.numeric(y_int) - out$fitted.values
  out$glm <- list(
    deviance      = g$deviance,
    null.deviance = g$null.deviance,
    df.null       = g$df.null,
    df.residual   = g$df.residual,
    aic           = g$aic,
    converged     = isTRUE(g$converged),
    iter          = g$iter,
    y.levels      = y_levels
  )
  out
}

# Internal: refit each bag replicate as binomial. We don't keep each bag's
# bootstrap rows — only its dirs/cuts/selected.terms. BUG-001 (v0.0.0.9029):
# the bootstrap indices are now captured in the *main* bagging loop and
# passed in via `idx_list`, so the GLM refit lands on the SAME rows the
# basis was selected on. Previously this function redrew from the live RNG
# and got a different bootstrap sample whenever `seed.cv = NULL`.
#
# @keywords internal
.ares_refit_boot_binomial <- function(fits, x_full, y_full, idx_list,
                                      weights = NULL) {
  n_full <- nrow(x_full)
  if (is.null(idx_list) || length(idx_list) != length(fits))
    stop("ares: internal error -- idx_list length mismatch in bag refit.")
  for (b in seq_along(fits)) {
    idx <- idx_list[[b]]
    x_b <- x_full[idx, , drop = FALSE]
    y_b <- y_full[idx]
    bx_b <- mars_basis_cpp(x_b, fits[[b]]$dirs, fits[[b]]$cuts,
                           as.integer(fits[[b]]$selected.terms))
    w_b <- if (!is.null(weights)) {
      ww <- weights[idx]; m <- mean(ww); if (m > 0) ww / m else ww
    } else rep(1, length(y_b))
    g <- tryCatch(
      suppressWarnings(stats::glm.fit(
        x = bx_b, y = y_b, weights = w_b,
        family = stats::binomial(), intercept = FALSE,
        control = list(maxit = 50L, epsilon = 1e-8, trace = FALSE))),
      error = function(e) NULL
    )
    if (!is.null(g)) {
      cf_b <- as.numeric(g$coefficients)
      if (any(!is.finite(cf_b))) cf_b[!is.finite(cf_b)] <- 0
      fits[[b]]$coefficients <- cf_b
    }
  }
  fits
}

# ============================================================================
#  Internal: post-hoc GLM refit (family = "poisson" or "gamma", log link)
# ============================================================================
#
# Same strategy as .ares_refit_binomial: forward + backward already picked
# terms via OLS; here we refit the chosen basis under poisson(link='log') or
# Gamma(link='log') via stats::glm.fit. The intercept column is included so
# we pass intercept = FALSE to glm.fit. Linear predictor is clamped to a
# wide band ([-50, 50] for poisson, [-20, 20] for gamma) so exp() stays
# finite under extrapolation; well-conditioned in-sample designs never hit
# the clamp.
#
# @keywords internal
.ares_refit_glm_loglink <- function(out, y, family_name, weights = NULL) {
  bx <- out$bx
  if (is.null(bx) || nrow(bx) == 0L)
    stop("ares: family='", family_name, "' but $bx is empty; cannot fit GLM.")
  fam_obj <- switch(family_name,
                    poisson = stats::poisson(link = "log"),
                    gamma   = stats::Gamma(link = "log"))
  # glm.fit starts from a family-specific default mu. For Gamma, default
  # starting mu = y can blow up the IRLS when y has order-of-magnitude
  # variation; passing an explicit mustart = y is fine here because y > 0
  # by validation upstream.
  must <- if (family_name == "gamma") y else NULL
  w_use <- if (is.null(weights)) rep(1, length(y)) else weights
  g <- suppressWarnings(stats::glm.fit(
    x = bx, y = y, family = fam_obj,
    weights = w_use, mustart = must,
    intercept = FALSE,
    control = list(maxit = 100L, epsilon = 1e-8, trace = FALSE)
  ))
  cf <- g$coefficients
  if (any(!is.finite(cf))) cf[!is.finite(cf)] <- 0
  names(cf) <- colnames(bx)
  lp <- drop(bx %*% cf)
  clamp_hi <- if (family_name == "poisson") 50 else 20
  clamp_lo <- if (family_name == "poisson") -50 else -20
  lp_clamped <- pmin(pmax(lp, clamp_lo), clamp_hi)
  out$coefficients     <- cf
  out$linear.predictor <- lp_clamped
  out$fitted.values    <- exp(lp_clamped)
  out$residuals        <- as.numeric(y) - out$fitted.values
  out$glm <- list(
    deviance      = g$deviance,
    null.deviance = g$null.deviance,
    df.null       = g$df.null,
    df.residual   = g$df.residual,
    aic           = g$aic,
    converged     = isTRUE(g$converged),
    iter          = g$iter,
    family        = family_name,
    link          = "log"
  )
  out
}

# Internal: refit each bag replicate under the log-link GLM family. Mirrors
# .ares_refit_boot_binomial. BUG-001 (v0.0.0.9029): bootstrap indices now
# come in via `idx_list` from the central bagging loop so the GLM refit
# uses the SAME rows the basis was selected on (previously, redrawing
# from the live RNG silently desynchronised under `seed.cv = NULL`).
#
# @keywords internal
.ares_refit_boot_glm_loglink <- function(fits, x_full, y_full, family_name,
                                          idx_list, weights = NULL) {
  n_full <- nrow(x_full)
  if (is.null(idx_list) || length(idx_list) != length(fits))
    stop("ares: internal error -- idx_list length mismatch in bag refit.")
  fam_obj <- switch(family_name,
                    poisson = stats::poisson(link = "log"),
                    gamma   = stats::Gamma(link = "log"))
  for (b in seq_along(fits)) {
    idx <- idx_list[[b]]
    x_b <- x_full[idx, , drop = FALSE]
    y_b <- y_full[idx]
    bx_b <- mars_basis_cpp(x_b, fits[[b]]$dirs, fits[[b]]$cuts,
                           as.integer(fits[[b]]$selected.terms))
    must <- if (family_name == "gamma") y_b else NULL
    w_b <- if (!is.null(weights)) {
      ww <- weights[idx]; m <- mean(ww); if (m > 0) ww / m else ww
    } else rep(1, length(y_b))
    g <- tryCatch(
      suppressWarnings(stats::glm.fit(
        x = bx_b, y = y_b, weights = w_b,
        family = fam_obj, intercept = FALSE,
        mustart = must,
        control = list(maxit = 100L, epsilon = 1e-8, trace = FALSE))),
      error = function(e) NULL
    )
    if (!is.null(g)) {
      cf_b <- as.numeric(g$coefficients)
      if (any(!is.finite(cf_b))) cf_b[!is.finite(cf_b)] <- 0
      fits[[b]]$coefficients <- cf_b
    }
  }
  fits
}

# ============================================================================
#  .ares_fit_varmod — fit the residual-variance model (gaussian only)
# ============================================================================
#
# Two flavours:
#   "const": store a single sigma_hat estimated from the (weighted) RSS
#            and an effective residual df. PIs are constant-width:
#               yhat +/- qt(1 - alpha/2, df) * sigma_hat.
#   "lm":    fit lm(|residual| ~ fitted) and multiply the prediction by
#            sqrt(pi/2) to convert E|N(0, sigma)| = sigma * sqrt(2/pi) back
#            to sigma. Captures simple heteroscedasticity; the model is
#            stored as a small list of intercept + slope (no need for the
#            full lm object).
#
# Inputs come from the fully-finalised gaussian fit. Both flavours store
# sigma estimates so predict.ares() can build PIs without re-running any
# costly step.
#
# @keywords internal
.ares_fit_varmod <- function(out, varmod, weights = NULL) {
  resid <- out$residuals
  yhat  <- out$fitted.values
  n  <- length(resid)
  k  <- length(out$coefficients)   # effective parameters used
  df <- max(n - k, 1L)             # residual df (guard against k >= n)
  if (varmod == "const") {
    if (is.null(weights)) {
      sigma_hat <- sqrt(sum(resid * resid) / df)
    } else {
      # Weighted residual variance: sum(w * r^2) / sum(w) is the natural
      # estimator under the WLS objective; df-adjust via (sum(w) - k).
      sw <- sum(weights)
      sigma_hat <- sqrt(sum(weights * resid * resid) / max(sw - k, 1))
    }
    list(type = "const", sigma_hat = sigma_hat, df = df)
  } else {                          # varmod == "lm"
    # |resid| ~ yhat (intercept + slope). lm() is overkill — closed form
    # gives the same numerical result for n moderate. Use weighted form
    # when weights are supplied.
    ar <- abs(resid)
    if (is.null(weights)) {
      m_y <- mean(yhat); m_r <- mean(ar)
      xc <- yhat - m_y;  rc <- ar  - m_r
      denom <- sum(xc * xc)
      slope <- if (denom > 0) sum(xc * rc) / denom else 0
      intercept <- m_r - slope * m_y
      sigma_const <- sqrt(pi / 2) * m_r
    } else {
      sw <- sum(weights)
      m_y <- sum(weights * yhat) / sw
      m_r <- sum(weights * ar)   / sw
      xc <- yhat - m_y; rc <- ar - m_r
      denom <- sum(weights * xc * xc)
      slope <- if (denom > 0) sum(weights * xc * rc) / denom else 0
      intercept <- m_r - slope * m_y
      sigma_const <- sqrt(pi / 2) * m_r
    }
    # E|N(0,sigma)| = sigma * sqrt(2/pi) -> sigma = |resid| * sqrt(pi/2).
    # BUG-003 (v0.0.0.9029) extrapolation floor: also stash an in-sample
    # range of yhat and a meaningful sigma fallback so predict() can warn
    # and floor (rather than silently collapse to ~1e-12 PI widths) when
    # the linear MAD model predicts a non-positive value at extrapolation
    # rows. `sigma_const_floor` is the mean-resid-derived constant sigma
    # (`varmod = "const"` equivalent on the same fit); we floor at
    # max(0.05 * sigma_const_floor, smallest positive in-sample MAD).
    yhat_min <- min(yhat)
    yhat_max <- max(yhat)
    # in-sample predicted MAD values (sigma scale)
    in_mad <- sqrt(pi / 2) * pmax(intercept + slope * yhat, 0)
    pos_in_mad <- in_mad[in_mad > 0]
    min_pos_in_mad <- if (length(pos_in_mad)) min(pos_in_mad) else sigma_const
    sigma_floor <- max(0.05 * sigma_const,
                       0.5  * min_pos_in_mad,
                       .Machine$double.eps)
    list(type = "lm",
         intercept = as.numeric(intercept),
         slope     = as.numeric(slope),
         scale     = sqrt(pi / 2),
         df        = df,
         yhat_min  = as.numeric(yhat_min),
         yhat_max  = as.numeric(yhat_max),
         sigma_const_floor = as.numeric(sigma_const),
         sigma_floor = as.numeric(sigma_floor))
  }
}

# ============================================================================
#  .term_labels — pretty-print MARS terms in earth's "h(x-cut)" format
# ============================================================================
#
# Walks the (dirs, cuts) matrices produced by the C++ engine and turns each
# row into a single human-readable label. Used for rownames of $bx, $dirs,
# $cuts, and as the names() of $coefficients.
#
# dirs[t, j] codes: 0 = variable unused, +/-1 = hinge sign, 2 = linear basis
# (auto.linpreds path). cuts[t, j] holds the knot for hinge terms; ignored
# when dirs[t, j] == 2.
#
# @keywords internal
.term_labels <- function(dirs, cuts, namesx) {
  M <- nrow(dirs)
  p <- ncol(dirs)
  labs <- character(M)
  labs[1] <- "(Intercept)"
  for (t in seq_len(M)[-1L]) {
    parts <- character(0)
    for (j in seq_len(p)) {
      d <- dirs[t, j]
      if (d == 0) next
      cv <- cuts[t, j]
      if (d == 2) {
        parts <- c(parts, namesx[j])
      } else if (d == 1) {
        parts <- c(parts, sprintf("h(%s-%g)", namesx[j], cv))
      } else if (d == -1) {
        parts <- c(parts, sprintf("h(%g-%s)", cv, namesx[j]))
      }
    }
    labs[t] <- if (length(parts)) paste(parts, collapse = " * ") else paste0("term", t)
  }
  labs
}

# ---- CV orchestration -------------------------------------------------------

# Internal: K-fold cross-validated subset-size selection for ares.
#
# Strategy:
#   1. For each repetition r in 1..ncross and each fold k in 1..nfold:
#      - Fit ares on the train rows with pmethod="backward" and return_path=1.
#        That returns the per-size backward path (each size's surviving terms
#        and OLS coefficients on the train fold).
#      - For each candidate size s, build the basis matrix on the *test rows*
#        using the fold fit's dirs/cuts/path.subsets[[s]], compute the
#        size-s prediction (bx_test %*% path.coefs[[s]]), and the size-s
#        holdout MSE.
#   2. Average MSE across folds and repetitions to get mean_MSE(s).
#   3. Pick size_star = argmin_s mean_MSE(s). Tie-break: smallest s.
#   4. Refit on full data with force_size = size_star.
#
# Notes:
#   - Each fold's forward pass is independent. A future optimisation could
#     reuse the fast.k cache across folds (see Phase 3 v0.20). v0.13 uses
#     plain repeated calls for clarity.
#   - Stratification: regression-only -- quantile-bin y into nfold bins and
#     round-robin within bins. ncross > 1 reseeds the partition.
# ============================================================================
#  .ares_cv_fit — K-fold CV pruning (pmethod = "cv")
# ============================================================================
#
# Drives the cross-validated subset-selection path:
#   1. Build stratified or random fold partitions (`nfold` folds repeated
#      `ncross` times).
#   2. For each (rep, fold), fit forward + backward on the training rows
#      using the C++ engine with `return_path = 1` so the backward path
#      (subset terms at every size 1..M) is returned.
#   3. Score every subset size on the held-out rows via `mars_basis_cpp` +
#      one matrix multiply per size — no extra fits.
#   4. Aggregate fold-MSE per size across the `nfold * ncross` evaluations.
#      Sizes not reached by every fold are excluded from the aggregate to
#      avoid the (smaller-MSE = more-terms) selection bias.
#   5. Pick `size_star` = argmin (or 1-SE rule if `cv.1se = TRUE`) and refit
#      the full data with `force_size = size_star`.
#
# Returns the full-data fit augmented with a $cv slot carrying:
#   nfold, ncross, stratify, cv.mse, cv.se, cv.n, cv.mse.mat, size, etc.
#
# @keywords internal
.ares_cv_fit <- function(x, y, degree, nk, penalty, thresh,
                         minspan, endspan, adjust_endspan, auto_linpreds,
                         fast_k, fast_beta, nprune, trace, nthreads,
                         nfold, ncross, stratify, seed_cv,
                         cv_1se = FALSE,
                         weights = NULL) {
  n <- nrow(x)
  if (n < nfold + 1L)
    stop("ares: nfold (", nfold, ") must be < n (", n, ").")

  # Stable, reproducible per-rep partitioning via a local RNG state save/
  # restore. We never advance the user's RNG when seed.cv is set.
  has_seed <- !is.null(seed_cv)
  if (has_seed) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(),
                       inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(as.integer(seed_cv))
  }

  make_folds <- function(rep_idx) {
    # Distinct partition per repetition. Quantile-stratified when stratify=TRUE.
    if (stratify) {
      # Quantile-bin y into ~nfold equal-count bins; within each bin shuffle
      # rows and assign fold IDs round-robin.
      ord <- order(y, sample(seq_len(n)))  # break ties at random
      fold_ids <- integer(n)
      fold_ids[ord] <- (seq_len(n) - 1L) %% nfold + 1L
    } else {
      perm <- sample(seq_len(n))
      fold_ids <- integer(n)
      fold_ids[perm] <- (seq_len(n) - 1L) %% nfold + 1L
    }
    fold_ids
  }

  # Run one fold and return per-size (weighted) MSE on holdout.
  # With weights supplied, the fold fit uses sample-weighted least squares
  # via mars_fit_cpp(weights_in = ...). Holdout scoring uses weighted-MSE
  # = sum(w_te * r^2) / sum(w_te) so model selection is consistent with the
  # WLS objective. When weights == NULL, both branches reduce to plain MSE.
  run_fold <- function(train_idx, test_idx) {
    x_tr <- x[train_idx, , drop = FALSE]
    y_tr <- y[train_idx]
    x_te <- x[test_idx,  , drop = FALSE]
    y_te <- y[test_idx]
    w_tr <- if (!is.null(weights)) {
      ww <- weights[train_idx]
      m <- mean(ww)
      if (m > 0) ww / m else ww
    } else NULL
    w_te <- if (!is.null(weights)) weights[test_idx] else NULL
    fit <- mars_fit_cpp(x_tr, as.numeric(y_tr), degree,
                        as.integer(nk %||% (min(200L, max(20L, 2L * ncol(x_tr))) + 1L)),
                        as.numeric(penalty %||% (if (degree > 1L) 3 else 2)),
                        thresh, minspan, endspan, adjust_endspan, auto_linpreds,
                        fast_k, fast_beta, as.integer(nprune %||% 0L),
                        0L, trace, nthreads, 0L, 1L,
                        weights_in = w_tr)
    M_full <- nrow(fit$dirs)
    out_mse <- rep(NA_real_, M_full)
    for (s in seq_len(M_full)) {
      subs <- fit$path.subsets[[s]]
      if (is.null(subs) || length(subs) == 0L) next
      coefs <- fit$path.coefs[[s]]
      if (length(coefs) != length(subs)) next  # defensive
      bx_te <- mars_basis_cpp(x_te, fit$dirs, fit$cuts,
                              as.integer(subs))
      yhat <- drop(bx_te %*% as.numeric(coefs))
      r <- y_te - yhat
      if (is.null(w_te)) {
        out_mse[s] <- mean(r * r)
      } else {
        sw <- sum(w_te)
        out_mse[s] <- if (sw > 0) sum(w_te * r * r) / sw else mean(r * r)
      }
    }
    out_mse
  }

  # Aggregate per-size MSE across nfold * ncross evaluations. We pad each
  # fold's vector to a common length (the maximum M observed) using NA, then
  # take the row-wise mean ignoring NAs.
  mse_list <- list()
  for (r in seq_len(ncross)) {
    fold_ids <- make_folds(r)
    for (k in seq_len(nfold)) {
      train_idx <- which(fold_ids != k)
      test_idx  <- which(fold_ids == k)
      if (length(test_idx) < 1L) next
      mse_list[[length(mse_list) + 1L]] <- run_fold(train_idx, test_idx)
    }
  }
  if (!length(mse_list))
    stop("ares: CV produced no evaluations -- check nfold / ncross.")
  max_M <- max(vapply(mse_list, length, integer(1L)))
  mse_mat <- matrix(NA_real_, nrow = length(mse_list), ncol = max_M)
  for (i in seq_along(mse_list)) {
    v <- mse_list[[i]]
    if (length(v) > 0L) mse_mat[i, seq_along(v)] <- v
  }
  cv_mse <- colMeans(mse_mat, na.rm = TRUE)
  cv_n   <- colSums(!is.na(mse_mat))
  full_n <- nrow(mse_mat)
  # Per-size CV-MSE standard error (across the nfold * ncross evaluations).
  # Only meaningful where cv_n > 1; treated as 0 otherwise.
  cv_se <- rep(NA_real_, length(cv_mse))
  for (s in seq_along(cv_mse)) {
    v <- mse_mat[, s]
    v <- v[is.finite(v)]
    if (length(v) > 1L) {
      cv_se[s] <- stats::sd(v) / sqrt(length(v))
    } else if (length(v) == 1L) {
      cv_se[s] <- 0
    }
  }
  # Prefer sizes achieved by ALL folds (otherwise mean is biased toward
  # folds whose forward pass happened to keep more terms). Fall back to
  # any-fold sizes only if no size was achieved by every fold.
  ok_all <- is.finite(cv_mse) & cv_n == full_n
  ok_any <- is.finite(cv_mse) & cv_n > 0
  if (any(ok_all)) {
    candidate_mask <- ok_all
  } else if (any(ok_any)) {
    candidate_mask <- ok_any
  } else {
    stop("ares: CV produced no usable size evaluations.")
  }
  cv_mse_finite <- cv_mse
  cv_mse_finite[!candidate_mask] <- Inf
  size_argmin <- which.min(cv_mse_finite)
  if (cv_1se) {
    # 1-SE rule: smallest size whose mean CV-MSE is within one standard
    # error of the minimum.
    se_at_min <- cv_se[size_argmin]
    if (!is.finite(se_at_min)) se_at_min <- 0
    threshold <- cv_mse_finite[size_argmin] + se_at_min
    in_thresh <- candidate_mask & is.finite(cv_mse) & cv_mse <= threshold
    if (any(in_thresh)) {
      size_star <- which(in_thresh)[1]   # smallest size meeting the rule
    } else {
      size_star <- size_argmin
    }
  } else {
    size_star <- size_argmin
  }

  # Final refit on full data with force_size = size_star. We still ask for
  # the path so we can return cv.mse to the caller.
  fit_full <- mars_fit_cpp(x, as.numeric(y), degree,
                           as.integer(nk %||% (min(200L, max(20L, 2L * ncol(x))) + 1L)),
                           as.numeric(penalty %||% (if (degree > 1L) 3 else 2)),
                           thresh, minspan, endspan, adjust_endspan,
                           auto_linpreds, fast_k, fast_beta,
                           as.integer(nprune %||% 0L),
                           0L, trace, nthreads,
                           as.integer(size_star), 1L,
                           weights_in = weights)
  # If the full-data forward pass produced fewer terms than size_star, the
  # C++ engine fell back to the GCV-best subset (force_size > M is ignored).
  # Detect and re-pick the best in-range size from cv_mse, then refit.
  M_full <- nrow(fit_full$dirs)
  if (size_star > M_full) {
    in_range <- candidate_mask & seq_along(cv_mse) <= M_full
    if (any(in_range)) {
      cv_mse_in <- cv_mse
      cv_mse_in[!in_range] <- Inf
      size_star <- which.min(cv_mse_in)
      fit_full <- mars_fit_cpp(x, as.numeric(y), degree,
                               as.integer(nk %||% (min(200L, max(20L, 2L * ncol(x))) + 1L)),
                               as.numeric(penalty %||% (if (degree > 1L) 3 else 2)),
                               thresh, minspan, endspan, adjust_endspan,
                               auto_linpreds, fast_k, fast_beta,
                               as.integer(nprune %||% 0L),
                               0L, trace, nthreads,
                               as.integer(size_star), 1L,
                               weights_in = weights)
    }
  }
  # Annotate.
  fit_full$cv <- list(
    nfold       = nfold,
    ncross      = ncross,
    stratify    = stratify,
    cv.1se      = cv_1se,
    cv.mse      = cv_mse,
    cv.se       = cv_se,
    cv.n        = cv_n,
    cv.mse.mat  = mse_mat,
    size.argmin = size_argmin,
    size.star   = size_star
  )
  fit_full
}

# Tiny helper: %||% (used only inside .ares_cv_fit for safe defaulting).
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Autotune (v0.15+) ------------------------------------------------------

# Internal: inner-CV grid search over (degree, penalty). For each cell, run
# K-fold CV with pmethod="backward" (GCV-based) per fold so that penalty
# directly shapes each fold's size pick. The cell's score is the mean
# holdout MSE across folds. Pick the smallest-score cell, then refit the
# winner on the full data using pmethod="backward" with that
# (degree, penalty). Return the winner augmented with $autotune carrying
# the grid summary.
#
# v0.15 grid:
#   degree  in {1, 2} (and 3 when nk_eff >= 21).
#   penalty in {0.5*d, 1.0*d, 2.0*d, 3.0*d} for each candidate degree d.
#
# Subsequent versions extend this grid (v0.16: nk; v0.17: autotune.speed;
# v0.18: n.boot). Helper signature stays stable.
# ============================================================================
#  .ares_autotune — hands-free hyperparameter tuning (autotune = TRUE)
# ============================================================================
#
# Picks (degree, penalty, nk, fast.k) by inner K-fold CV grid search:
#
#   - Grid construction: degree in {1, 2} (plus 3 when nk allows);
#     penalty in {0.5, 1.0, 2.0, 3.0} * degree; nk in {1, 2, 4} * nk_default
#     (capped at 200; capped at 2x on high-p where nk_default >= 31);
#     fast.k driven by `autotune.speed`. Cells share one fold partition.
#
#   - Successive halving: after fold 1, any cell whose fold-1 MSE exceeds
#     1.5 * running best is eliminated; remaining folds are skipped for it.
#
#   - Warm-start (autotune.warmstart = TRUE): a 15% subsample (capped at
#     200 rows, gated by n >= 200) runs the grid first. If one degree's
#     best cell beats the next-best degree's by >5% / >20%, that cell is
#     used directly to refit on the full data — skipping the full-grid CV
#     entirely. Cuts wall-clock by ~5x on well-separated DGPs.
#
#   - Shared forward pass (Phase 3): cells grouped by (degree, nk, fast.k)
#     run one C++ forward per group per fold; backward replay via
#     `mars_backward_only_cpp` produces the per-cell subset paths.
#
# Returns a full-data ares fit with $autotune containing the grid table,
# the winning row, and the chosen hyperparameters.
#
# @keywords internal
.ares_autotune <- function(x, y, nk, thresh, minspan, endspan,
                           adjust_endspan, auto_linpreds, fast_k, fast_beta,
                           nprune, nfold, ncross, stratify, seed_cv,
                           cv_1se, trace, nthreads,
                           autotune_speed = "balanced",
                           warmstart = TRUE,
                           weights = NULL) {
  n <- nrow(x)
  p <- ncol(x)

  # Effective nk used for the grid (unless user pinned it).
  nk_eff <- if (length(nk) == 0L || is.na(nk[1])) {
    min(200L, max(20L, 2L * p)) + 1L
  } else as.integer(nk)

  # Build degree grid: include 3 only if nk_eff >= 31 (need plenty of forward
  # slots for a 3-way hinge to actually appear in the forward pass; with
  # nk = 2*p + 1 ~ 21 there's barely room past depth-2 selection).
  deg_grid <- c(1L, 2L)
  if (nk_eff >= 31L) deg_grid <- c(deg_grid, 3L)

  # (degree, penalty, nk) cells. nk multipliers are 1, 2, 4 of the default,
  # capped at 200. The default value (nk_default) is the user's nk if pinned,
  # else min(200, max(20, 2*p)) + 1.
  # ---- Warm-start (v0.19+) ----------------------------------------------
  # Pre-fit autotune on a 15% subsample (capped at 200 rows). If the
  # subsample's best-per-degree CV-MSE is at least 5% below the next-best
  # degree's best cell, adopt that degree's winner directly and skip
  # the full-data grid. Disabled when n < 200 (subsample would be tiny).
  warm_winner <- NULL
  if (warmstart && n >= 200L) {
    n_sub <- min(200L, max(100L, as.integer(round(0.20 * n))))
    has_seed <- !is.null(seed_cv)
    if (has_seed) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed_ws <- get(".Random.seed", envir = globalenv(),
                           inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed_ws, envir = globalenv()),
                add = TRUE)
      } else {
        on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
      }
      set.seed(as.integer(seed_cv) + 503L)
    }
    sub_idx <- sample.int(n, n_sub, replace = FALSE)
    x_sub <- x[sub_idx, , drop = FALSE]
    y_sub <- y[sub_idx]
    w_sub <- if (!is.null(weights)) {
      ww <- weights[sub_idx]
      m <- mean(ww)
      if (m > 0) ww / m else ww
    } else NULL
    sub_fit <- tryCatch(
      .ares_autotune(
        x = x_sub, y = y_sub,
        nk = as.integer(nk_eff), thresh = thresh,
        minspan = minspan, endspan = endspan,
        adjust_endspan = adjust_endspan, auto_linpreds = auto_linpreds,
        fast_k = fast_k, fast_beta = fast_beta,
        nprune = nprune,
        nfold = 3L,             # cheap inner CV
        ncross = 1L,
        stratify = stratify,
        seed_cv = if (has_seed) as.integer(seed_cv) + 503L else NULL,
        cv_1se = cv_1se, trace = trace, nthreads = nthreads,
        autotune_speed = "fast",   # fast.k = 5 always inside subsample
        warmstart = FALSE,
        weights = w_sub
      ),
      error = function(e) NULL
    )
    if (!is.null(sub_fit) && !is.null(sub_fit$autotune)) {
      g <- sub_fit$autotune$grid
      g_alive <- g[is.finite(g$cv_mse), , drop = FALSE]
      if (nrow(g_alive) >= 1L) {
        best_per_deg <- tapply(g_alive$cv_mse, g_alive$degree,
                                min, na.rm = TRUE)
        decisive <- FALSE
        if (length(best_per_deg) == 1L) {
          decisive <- TRUE
        } else {
          ord <- order(best_per_deg)
          best_score   <- as.numeric(best_per_deg[ord[1]])
          runner_score <- as.numeric(best_per_deg[ord[2]])
          # Require >=10% gap on subsample to commit to early-exit; 5% was
          # empirically too lax (committed to deg=1 on a 75-row subsample of
          # a friedman-1 design where deg=2 wins on full data).
          if (is.finite(best_score) && is.finite(runner_score) &&
              best_score < 0.90 * runner_score) decisive <- TRUE
        }
        if (decisive) {
          warm_winner <- list(
            degree   = sub_fit$autotune$degree,
            penalty  = sub_fit$autotune$penalty,
            nk       = sub_fit$autotune$nk,
            fast_k   = sub_fit$autotune$fast_k
          )
        }
      }
    }
  }
  if (!is.null(warm_winner)) {
    fit_full <- mars_fit_cpp(
      x, y, as.integer(warm_winner$degree),
      as.integer(warm_winner$nk), as.numeric(warm_winner$penalty),
      thresh, minspan, endspan, adjust_endspan, auto_linpreds,
      as.integer(warm_winner$fast_k), fast_beta,
      if (length(nprune) == 0L || is.na(nprune[1])) {
        as.integer(warm_winner$nk)
      } else as.integer(nprune),
      0L, trace, nthreads, 0L, 0L,
      weights_in = weights
    )
    fit_full$autotune <- list(
      grid          = data.frame(
        degree     = integer(0), penalty = numeric(0),
        nk         = integer(0), fast_k  = integer(0),
        cv_mse     = numeric(0), eliminated = logical(0)
      ),
      best          = NA_integer_,
      degree        = warm_winner$degree,
      penalty       = warm_winner$penalty,
      nk            = warm_winner$nk,
      fast_k        = warm_winner$fast_k,
      cv_mse        = NA_real_,
      nfold         = nfold,
      ncross        = ncross,
      stratify      = stratify,
      n_eliminated  = 0L,
      speed         = autotune_speed,
      warmstart     = TRUE
    )
    return(fit_full)
  }

  # v0.0.0.9021: Cap nk_grid at 2x of nk_eff when nk_eff is already big
  # (>=31, i.e. p>=15). The 4x multiplier (=>nk approaching 200 for p=20)
  # produces forwards that emit ~M=130+ terms; backward-replay cost on the
  # autotune CV grid scales as O(M^4) per cell, so capping M downstream
  # gives ~5x wall-clock savings without an empirical MSE penalty on the
  # mlbench highdim sweep. For low-p problems (nk_eff < 31, p < 15) keep
  # the full c(1x, 2x, 4x) sweep -- the cells are cheap and the extra
  # multiplier occasionally pays in MSE.
  nk_grid <- if (nk_eff >= 31L) {
    unique(pmin(c(nk_eff, 2L * nk_eff), 200L))
  } else {
    unique(pmin(c(nk_eff, 2L * nk_eff, 4L * nk_eff), 200L))
  }

  # autotune_speed determines fast.k policy:
  #   "quality"  : fast.k = 0 always (no cache, slowest, most accurate).
  #   "fast"     : fast.k = 5 always (aggressive cache, cheapest).
  #   "balanced" : sweep fast.k in {10, 25} (and 0, == "no cache", on
  #                low-p only) inside the grid; pick smallest within 1%
  #                of best.
  #
  # v0.0.0.9022: on high-p (nk_eff >= 31), drop fast.k = 0 from the balanced
  # sweep. Cache-disabled forwards at deg>=2 cost ~3-5x the cached ones, and
  # empirically the within-1% set on highdim DGPs never contains a fast.k=0
  # cell -- the cache prevents over-detailed knot scans that slightly hurt
  # CV-MSE. Net: 1.5x extra wall-clock cut on top of the nk-grid cap, with
  # MSE preserved or improved.
  fk_grid <- switch(autotune_speed,
                    quality  = 0L,
                    fast     = 5L,
                    balanced = if (nk_eff >= 31L) c(10L, 25L) else c(10L, 25L, 0L))

  cells <- list()
  for (d in deg_grid)
    for (mult in c(0.5, 1.0, 2.0, 3.0))
      for (nk_c in nk_grid)
        for (fk in fk_grid)
          cells[[length(cells) + 1L]] <- list(degree  = d,
                                              penalty = mult * d,
                                              nk      = as.integer(nk_c),
                                              fast_k  = as.integer(fk))

  # Build a single shared fold partition (function of seed.cv + stratify) so
  # all cells score on the same folds. This is critical for fair comparison
  # between cells. Stratification is regression-only (quantile-binned y).
  has_seed <- !is.null(seed_cv)
  if (has_seed) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(),
                       inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(as.integer(seed_cv))
  }
  fold_lists <- vector("list", ncross)
  for (r in seq_len(ncross)) {
    if (stratify) {
      ord <- order(y, sample(seq_len(n)))
      fid <- integer(n)
      fid[ord] <- (seq_len(n) - 1L) %% nfold + 1L
    } else {
      perm <- sample(seq_len(n))
      fid <- integer(n)
      fid[perm] <- (seq_len(n) - 1L) %% nfold + 1L
    }
    fold_lists[[r]] <- fid
  }

  # Score one fold for one cell. Returns NA on failure.
  # Currently unused by the shared-forward fast path; kept for clarity/test
  # comparison. Weights composed identically to the shared-forward path.
  score_one <- function(d, pen, nk_c, fk, tr, te) {
    w_tr <- if (!is.null(weights)) {
      ww <- weights[tr]; m <- mean(ww); if (m > 0) ww / m else ww
    } else NULL
    w_te <- if (!is.null(weights)) weights[te] else NULL
    fit_k <- mars_fit_cpp(
      x[tr, , drop = FALSE], as.numeric(y[tr]),
      as.integer(d), as.integer(nk_c), as.numeric(pen),
      thresh, minspan, endspan, adjust_endspan, auto_linpreds,
      as.integer(fk), fast_beta,
      if (length(nprune) == 0L || is.na(nprune[1])) as.integer(nk_c) else as.integer(nprune),
      0L, trace, nthreads, 0L, 0L,
      weights_in = w_tr
    )
    bx_te <- mars_basis_cpp(x[te, , drop = FALSE],
                            fit_k$dirs, fit_k$cuts,
                            as.integer(fit_k$selected.terms))
    yhat <- drop(bx_te %*% as.numeric(fit_k$coefficients))
    r_te <- y[te] - yhat
    if (is.null(w_te)) mean(r_te * r_te)
    else {
      sw <- sum(w_te); if (sw > 0) sum(w_te * r_te * r_te) / sw
      else mean(r_te * r_te)
    }
  }

  # Successive-halving + shared forward pass (v0.20+).
  # Cells are grouped by (degree, nk, fast.k): the forward pass output is
  # independent of penalty, so for each fold we run forward ONCE per group
  # and call mars_backward_only_cpp() per cell within the group with its
  # specific penalty. mars_backward_only_cpp builds B from dirs/cuts on the
  # fold's training rows and replays GCV-backward -- same numerical path as
  # mars_fit_cpp's backward block. Halving factor 1.5 (eliminate cells whose
  # fold-1 MSE exceeds 1.5 x running_best).
  halving_factor <- 1.5
  ncell <- length(cells)
  fold_pairs <- list()
  for (r in seq_len(ncross))
    for (k in seq_len(nfold))
      fold_pairs[[length(fold_pairs) + 1L]] <- c(r, k)

  # Group cells by (degree, nk, fast_k). cell_group[i] gives the group index.
  group_key <- vapply(cells, \(c)
    sprintf("%d|%d|%d", c$degree, c$nk, c$fast_k), character(1L))
  group_levels <- unique(group_key)
  cell_group <- match(group_key, group_levels)
  groups <- vector("list", length(group_levels))
  for (g in seq_along(group_levels)) {
    members <- which(cell_group == g)
    first <- cells[[members[1]]]
    groups[[g]] <- list(degree = first$degree, nk = first$nk,
                         fast_k = first$fast_k,
                         members = members)
  }

  fold_mse <- matrix(NA_real_, nrow = ncell, ncol = length(fold_pairs))
  alive <- rep(TRUE, ncell)
  for (fp_i in seq_along(fold_pairs)) {
    r_idx <- fold_pairs[[fp_i]][1]
    k_idx <- fold_pairs[[fp_i]][2]
    fid <- fold_lists[[r_idx]]
    tr <- which(fid != k_idx); te <- which(fid == k_idx)
    if (length(te) < 1L) next
    x_tr <- x[tr, , drop = FALSE]
    y_tr <- as.numeric(y[tr])
    x_te <- x[te, , drop = FALSE]
    y_te <- y[te]
    # Fold-level weights: rescale train weights to mean 1 (matches the C++
    # engine convention); test weights left raw (we score weighted MSE).
    w_tr <- if (!is.null(weights)) {
      ww <- weights[tr]; m <- mean(ww); if (m > 0) ww / m else ww
    } else NULL
    w_te <- if (!is.null(weights)) weights[te] else NULL

    for (g in seq_along(groups)) {
      gi <- groups[[g]]
      # Skip group if all members are eliminated.
      if (!any(alive[gi$members])) next

      # Use the smallest penalty (within members) for the forward fit so the
      # forward pass terminates only on rss-improvement / nk, not on GCV
      # (which doesn't gate the forward anyway). Penalty has no effect on
      # forward, so any value works.
      pen0 <- min(vapply(cells[gi$members], \(c) c$penalty, numeric(1L)))
      ff <- tryCatch(
        mars_fit_cpp(
          x_tr, y_tr, as.integer(gi$degree), as.integer(gi$nk),
          as.numeric(pen0), thresh, minspan, endspan,
          adjust_endspan, auto_linpreds,
          as.integer(gi$fast_k), fast_beta,
          if (length(nprune) == 0L || is.na(nprune[1])) {
            as.integer(gi$nk)
          } else as.integer(nprune),
          0L,           # pmethod=backward (need a forward pass)
          trace, nthreads, 0L, 0L,
          weights_in = w_tr
        ),
        error = function(e) NULL
      )
      if (is.null(ff)) {
        for (mi in gi$members)
          if (alive[mi]) fold_mse[mi, fp_i] <- NA_real_
        next
      }
      dirs_full <- ff$dirs
      cuts_full <- ff$cuts

      # Per-cell backward replay with its own penalty. Same fold weights
      # are passed through so the IRLS/WLS objective stays consistent.
      for (mi in gi$members) {
        if (!alive[mi]) next
        ce <- cells[[mi]]
        bk <- tryCatch(
          mars_backward_only_cpp(
            x_tr, y_tr, dirs_full, cuts_full,
            as.numeric(ce$penalty),
            if (length(nprune) == 0L || is.na(nprune[1])) {
              as.integer(ce$nk)
            } else as.integer(nprune),
            nthreads, 0L, 0L,
            weights_in = w_tr
          ),
          error = function(e) NULL
        )
        if (is.null(bk)) { fold_mse[mi, fp_i] <- NA_real_; next }
        bx_te <- mars_basis_cpp(x_te, dirs_full, cuts_full,
                                as.integer(bk$selected.terms))
        yhat <- drop(bx_te %*% as.numeric(bk$coefficients))
        r_te <- y_te - yhat
        if (is.null(w_te)) {
          fold_mse[mi, fp_i] <- mean(r_te * r_te)
        } else {
          sw <- sum(w_te)
          fold_mse[mi, fp_i] <- if (sw > 0) sum(w_te * r_te * r_te) / sw
                                 else mean(r_te * r_te)
        }
      }
    }

    # Halving check after fold 1 only.
    if (fp_i == 1L) {
      f1 <- fold_mse[, 1]
      ok <- alive & is.finite(f1)
      if (any(ok)) {
        run_best <- min(f1[ok])
        cutoff <- halving_factor * run_best
        elim <- alive & (!is.finite(f1) | f1 > cutoff)
        if (any(elim)) alive[elim] <- FALSE
      }
    }
  }

  # Cell score: mean over the folds that actually got run (alive cells get
  # all folds; eliminated cells get only fold 1).
  scores <- numeric(ncell)
  for (i in seq_len(ncell)) {
    v <- fold_mse[i, ]
    v <- v[is.finite(v)]
    scores[i] <- if (length(v)) mean(v) else Inf
  }
  if (!any(is.finite(scores)))
    stop("ares: autotune found no finite CV-MSE on any grid cell.")

  # Tie-break: smallest score; on tie, prefer smaller degree (parsimony),
  # then smaller penalty.
  # Cell order key for default (non-balanced) tie-break.
  ord_cells <- order(scores,
                     vapply(cells, \(c) c$degree, integer(1L)),
                     vapply(cells, \(c) c$nk, integer(1L)),
                     vapply(cells, \(c) c$penalty, numeric(1L)),
                     vapply(cells, \(c) c$fast_k, integer(1L)))

  # Balanced mode: among cells within 1% of the argmin score, prefer the
  # smallest non-zero fast_k (= fastest cache setting). fast_k = 0 (cache
  # disabled) gets picked only if it is strictly the global best within
  # 1% (i.e. nothing else qualifies).
  if (autotune_speed == "balanced" && any(is.finite(scores))) {
    s_min <- min(scores[is.finite(scores)])
    within1 <- is.finite(scores) & scores <= s_min * 1.01
    if (any(within1)) {
      # Sort within1 cells by: fastest fast_k first (smaller non-zero fast_k
      # is cheaper). Treat fast_k = 0 as the highest sort key (only chosen
      # when no positive-fast_k cell qualifies).
      fk_vec <- vapply(cells, \(c) c$fast_k, integer(1L))
      score_key <- ifelse(fk_vec == 0L, .Machine$integer.max, fk_vec)
      candidates <- which(within1)
      ord_within <- candidates[order(score_key[candidates],
                                     scores[candidates],
                                     vapply(cells[candidates],
                                            \(c) c$degree,
                                            integer(1L)),
                                     vapply(cells[candidates],
                                            \(c) c$nk,
                                            integer(1L)))]
      # Override the argmin pick with the balanced-rule winner.
      ord_cells <- c(ord_within, setdiff(ord_cells, ord_within))
    }
  }
  best <- ord_cells[1]
  best_d   <- cells[[best]]$degree
  best_pen <- cells[[best]]$penalty

  best_nk  <- cells[[best]]$nk
  best_fk  <- cells[[best]]$fast_k
  # Refit winner on the full data, GCV-backward, with the chosen
  # (degree, penalty, nk, fast.k).
  fit_full <- mars_fit_cpp(
    x, y, as.integer(best_d), as.integer(best_nk), as.numeric(best_pen),
    thresh, minspan, endspan, adjust_endspan, auto_linpreds,
    as.integer(best_fk), fast_beta,
    if (length(nprune) == 0L || is.na(nprune[1])) as.integer(best_nk) else as.integer(nprune),
    0L, trace, nthreads, 0L, 0L,
    weights_in = weights
  )

  grid_df <- data.frame(
    degree     = vapply(cells, \(c) c$degree, integer(1L)),
    penalty    = vapply(cells, \(c) c$penalty, numeric(1L)),
    nk         = vapply(cells, \(c) c$nk, integer(1L)),
    fast_k     = vapply(cells, \(c) c$fast_k, integer(1L)),
    cv_mse     = scores,
    eliminated = !alive,
    stringsAsFactors = FALSE
  )
  fit_full$autotune <- list(
    grid          = grid_df,
    best          = best,
    degree        = best_d,
    penalty       = best_pen,
    nk            = best_nk,
    fast_k        = best_fk,
    cv_mse        = scores[best],
    nfold         = nfold,
    ncross        = ncross,
    stratify      = stratify,
    n_eliminated  = sum(!alive),
    speed         = autotune_speed,
    warmstart     = FALSE
  )
  fit_full
}

