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
#' Fits a MARS model (Friedman 1991) using a fast least-squares forward pass
#' with a parallelized knot search and a backward subset-selection that
#' minimizes the GCV criterion. The implementation aims for numerical parity
#' with [earth::earth()] on the gaussian-only core while taking advantage of
#' multi-core CPUs to reduce wall-clock fitting time.
#'
#' Two interfaces are provided. `ares.default()` accepts a numeric predictor
#' matrix `x` and a numeric response vector `y`. `ares.formula()` accepts a
#' model formula and a data frame.
#'
#' @param x A numeric matrix or data frame of predictors, OR a model formula
#'   when calling the formula method.
#' @param y A numeric response vector. Length must equal `nrow(x)`.
#' @param data Used by the formula method only. A data frame.
#' @param degree Maximum interaction degree. Default 1.
#' @param nk Maximum number of basis terms in the forward pass. Default
#'   `min(200, max(20, 2 * ncol(x))) + 1`.
#' @param penalty GCV penalty. Default `if (degree > 1) 3 else 2`.
#' @param thresh Forward-pass relative-RSS early-stop threshold. Default 0.001.
#' @param minspan Minimum knot span. 0 (default) selects an automatic value.
#' @param endspan Knot offset from the data ends. 0 (default) selects automatic.
#' @param adjust.endspan Multiplier applied to `endspan` when the candidate
#'   hinge would deepen an existing interaction (parent term has
#'   degree at least 1). Default `1` (no adjustment). Earth's analogous
#'   default is `2`; ares's heuristic is empirically mixed on the
#'   inst/sims grid (helps a few cells, hurts others), so the default
#'   is conservative. Pass `2` to enable.
#' @param auto.linpreds If `TRUE`, and the best forward-pass knot for a
#'   candidate hinge sits at the boundary of its variable's eligible range,
#'   substitutes a linear term (`dirs = 2`) for the hinge pair. **Default
#'   `FALSE`.** Earth's analogous default is `TRUE`, but earth applies a
#'   stricter linearity test that ares's "boundary knot" heuristic only
#'   approximates; current ares implementation regresses parity on most
#'   benchmark cells. Treat as experimental.
#' @param fast.k Fast-MARS priority cache size (Friedman 1993). The forward
#'   pass always rescores pairs whose parent was added in the previous step;
#'   the top `fast.k` of the remaining stale pairs (ranked by age-discounted
#'   cached score) are also rescored, and the rest contribute their cached
#'   reduction. `0` disables the cache (every pair is rescored every step).
#'   Default `10`.
#' @param fast.beta Age-penalty for the fast-MARS priority cache: a stale
#'   pair's effective score is `cached_rss_red / (1 + fast.beta * age)`.
#'   Default `1.0`. Only matters when `fast.k > 0`.
#' @param nprune Maximum number of terms after backward pruning (default `nk`).
#' @param pmethod Pruning method: `"backward"` (default), `"none"`, or
#'   `"cv"` for K-fold cross-validated subset selection. `"cv"` requires
#'   `nfold > 0` (or implicitly sets `nfold = 5` when not specified).
#' @param nfold Number of CV folds (default `0` = no CV; GCV-based backward
#'   selection). When `nfold > 0`, `pmethod` is promoted to `"cv"`.
#' @param ncross Number of CV repetitions when `nfold > 0`. Each repetition
#'   builds a fresh fold partition; per-size mean MSE is averaged across
#'   `nfold * ncross` evaluations. Default `1`.
#' @param stratify If `TRUE` (default), folds are quantile-stratified on `y`
#'   so each fold has a similar response distribution. Useful on small `n` /
#'   skewed responses; cheap.
#' @param seed.cv Optional integer seed for fold partitioning. `NULL`
#'   (default) uses the current RNG state; pass an integer for reproducible
#'   CV partitions.
#' @param cv.1se If `TRUE`, applies the one-standard-error rule for
#'   `pmethod="cv"`: among sizes whose mean CV-MSE is within one
#'   standard error of the minimum, pick the smallest. Trades a bit of
#'   accuracy for parsimony. Default `FALSE` (pick `argmin` of mean
#'   CV-MSE, i.e. the size that minimises holdout error directly).
#' @param autotune If `TRUE`, runs an inner cross-validated grid search
#'   over `(degree, penalty)` and refits the winner. The grid is
#'   `degree in {1, 2}` (extended to `{1, 2, 3}` when `nk` is large
#'   enough) crossed with `penalty in {0.5*d, 1.0*d, 2.0*d, 3.0*d}` at
#'   each degree `d`. Inner CV uses `nfold` (default 5 when autotune
#'   is on and `nfold == 0`) and respects `seed.cv`. The chosen
#'   `(degree, penalty)` and the full grid CV-MSE are returned in
#'   `$autotune`. Default `FALSE` (current grf-style fast default
#'   path). Subsequent versions extend the grid (v0.16: `nk`,
#'   v0.17: `autotune.speed`, v0.18: `n.boot`).
#' @param autotune.speed One of `"balanced"` (default), `"quality"`,
#'   or `"fast"`. Only meaningful when `autotune = TRUE`.
#'   - `"quality"`: forces `fast.k = 0` for every cell (no fast-MARS
#'     priority cache; every (parent, var) pair rescored every step).
#'     Matches v0.10 quality at v0.10 wall-clock cost.
#'   - `"fast"`: forces `fast.k = 5`. Most aggressive cache; cheapest.
#'   - `"balanced"` (default): sweeps `fast.k in {10, 25, Inf}` (where
#'     `Inf` means "no cache") inside the autotune grid and picks the
#'     smallest `fast.k` whose mean CV-MSE is within 1% of the best.
#'     Trades a marginal accuracy hit for a marginal speed gain.
#' @param autotune.warmstart If `TRUE` (default; only meaningful when
#'   `autotune = TRUE`) and `n >= 200`, ares first runs autotune on a
#'   15% subsample (capped at 200 rows). If the subsample has a
#'   *decisive* winner -- best-per-degree CV-MSE more than 5% below the
#'   next-best degree's best cell -- that `(degree, penalty, nk, fast.k)`
#'   is used directly to refit on the full data, skipping the full-data
#'   grid. Cuts autotune wall-clock by ~5x on well-separated DGPs.
#' @param n.boot Number of bootstrap replicate fits for bagging. Default
#'   `0` (no bagging). When `> 0`, the result holds a list of `n.boot` ares
#'   fits in `$boot$fits`, each fit on a row-bootstrap sample of the data.
#'   The augmented `predict()` averages predictions across replicates plus
#'   the central fit and reports per-prediction standard deviation in
#'   `attr(predictions, "sd")`. Bagging composes with `autotune` --
#'   each replicate fits with the central fit's chosen
#'   `(degree, penalty, nk, fast.k)` so the cost is `(n.boot + 1) * fit_cost`,
#'   not the autotune grid times `n.boot`.
#' @param na.action Strategy for missing values in `x`. Either `"impute"`
#'   (default) or `"omit"`. `"impute"` replaces each column's `NA`s with
#'   that column's median (computed from the non-missing rows) and stores
#'   those medians on the fit so `predict()` can reapply them to future
#'   newdata. `"omit"` drops every row that has any `NA` in `x`. Both
#'   actions emit a warning summarising the affected rows / columns.
#'   Missing values in `y` are always a hard error -- there is no sensible
#'   way to impute the regression target. `NaN` and `+/-Inf` in `x` are
#'   also rejected outright regardless of `na.action`.
#' @param family Response family. Either `"gaussian"` (default; numeric `y`,
#'   identity link, OLS coefficients) or `"binomial"` (binary `y`, logit
#'   link). When `family = "binomial"`, `y` must be either a 0/1 numeric
#'   vector or a 2-level factor (any other shape is a hard error). The
#'   forward and backward passes still run on the original numeric `y` as
#'   if gaussian (matching earth's GLM strategy); after backward pruning
#'   the selected basis is refit on the binary response via
#'   `stats::glm.fit()` with `family = binomial()`, and `$coefficients`,
#'   `$fitted.values`, and `$linear.predictor` come from that GLM. Bagging,
#'   CV pruning, and autotune compose with `family = "binomial"`; the inner
#'   model-selection still uses Gaussian MSE on the latent scale (cheap and
#'   adequate for term selection — only the final coefficients use IRLS).
#' @param trace Trace level. 0 = silent (default), 1 = forward-pass progress.
#' @param nthreads Number of threads. 0 (default) selects
#'   `RcppParallel::defaultNumThreads()`. CRAN-distributed examples and the
#'   bundled vignette cap this at 2.
#' @param ... Additional arguments. Currently ignored.
#' @return An object of class `"ares"` -- a list with components
#'   `coefficients`, `bx`, `dirs`, `cuts`, `selected.terms`, `rss`, `gcv`,
#'   `rss.per.subset`, `gcv.per.subset`, `fitted.values`, `residuals`, `namesx`,
#'   `call`, plus echoed control parameters. When `family = "binomial"` the
#'   fit additionally carries `$family = "binomial"`, `$glm` (a small list
#'   with `deviance`, `null.deviance`, `df.null`, `df.residual`, `aic`,
#'   `converged`, `iter`, `y.levels`), `$linear.predictor` (length `n`),
#'   and `$fitted.values` holds response-scale probabilities (`plogis(lp)`).
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics* 19(1):1-67.
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
                         family = c("gaussian", "binomial"),
                         trace = 0L, nthreads = 0L, ...) {
  cl <- match.call()
  pmethod <- match.arg(pmethod)
  autotune_speed <- match.arg(autotune.speed)
  na.action <- match.arg(na.action)
  family <- match.arg(family)
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
  if (is.data.frame(x)) {
    is_cat <- vapply(x, \(z) is.factor(z) || is.character(z), logical(1L))
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
        expanded_names = colnames(xmm)      # expanded dummy column names
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
      warmstart = isTRUE(autotune.warmstart)
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
                        cv_1se = isTRUE(cv.1se))
  } else {
    out <- mars_fit_cpp(x, as.numeric(y), degree, nk, penalty, thresh,
                        minspan, endspan, adjust_endspan, auto_linpreds,
                        fast_k, fast_beta,
                        nprune, pmethod_int, trace, nthreads_eff,
                        force_size = 0L, return_path = 0L)
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
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n_full, n_full, replace = TRUE)
      x_b <- x[idx, , drop = FALSE]
      y_b <- as.numeric(y)[idx]
      f_b <- mars_fit_cpp(
        x_b, y_b, as.integer(deg_b), as.integer(nk_b),
        as.numeric(pen_b),
        thresh, minspan, endspan, adjust_endspan, auto_linpreds,
        as.integer(fk_b), fast_beta,
        if (length(nprune) == 0L || is.na(nprune[1])) as.integer(nk_b) else as.integer(nprune),
        pmethod_int_b, trace, nthreads_eff, 0L, 0L
      )
      # Strip heavy fields that bagging doesn't need (bx is n_full x M).
      f_b$bx <- NULL
      boot_fits[[b]] <- f_b
    }
    out$boot <- list(
      fits   = boot_fits,
      n.boot = n_boot
    )
  }

  # Persist the NA-handling metadata so predict.ares() can reapply the
  # same imputation to newdata. `na.medians` is NULL when the training
  # data had no NAs or `na.action = "omit"` was used.
  out$na.action <- na.action
  out$na.medians <- na_medians
  # Factor-expansion metadata: NULL when training x was already a numeric
  # matrix (or a data.frame with no factor/character columns). When set,
  # predict.ares() replays the same model.matrix expansion on newdata.
  out$factor_info <- factor_info

  # ---- family bookkeeping --------------------------------------------------
  out$family <- family
  if (family == "binomial") {
    # Post-hoc GLM refit on the selected basis (earth strategy). The
    # forward + backward already picked terms via OLS on the binary y; now
    # IRLS estimates logit-scale coefficients on the same design. We use
    # `glm.fit` (no formula, no model.frame) for speed.
    out <- .ares_refit_binomial(out, y_int = as.integer(y),
                                y_levels = y_levels)
    # Refit each bag replicate too (each replicate already holds its own
    # selected basis on its bootstrap rows; we re-evaluate on those rows
    # and IRLS-refit). The original x/y are needed because bag fits only
    # carry dirs/cuts/selected.terms (not their bootstrap data).
    if (!is.null(out$boot) && length(out$boot$fits)) {
      out$boot$fits <- .ares_refit_boot_binomial(out$boot$fits,
                                                  x_full = x,
                                                  y_full = as.integer(y),
                                                  seed_cv = seed.cv)
    }
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
.ares_refit_binomial <- function(out, y_int, y_levels) {
  bx <- out$bx
  if (is.null(bx) || nrow(bx) == 0L)
    stop("ares: family='binomial' but $bx is empty; cannot fit GLM.")
  # The intercept column is the first basis term (all ones), so glm.fit gets
  # the full design and we pass `intercept = FALSE` to avoid duplicating it.
  g <- stats::glm.fit(x = bx, y = y_int, family = stats::binomial(),
                      intercept = FALSE)
  cf <- g$coefficients
  names(cf) <- colnames(bx)
  out$coefficients   <- cf
  out$linear.predictor <- drop(bx %*% cf)
  out$fitted.values  <- stats::plogis(out$linear.predictor)
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
# bootstrap rows — only its dirs/cuts/selected.terms. So we rebuild the bag's
# basis on a fresh bootstrap sample (using the same seed offset as the
# original bagging loop) and IRLS-fit it against that sample's binary y.
#
# @keywords internal
.ares_refit_boot_binomial <- function(fits, x_full, y_full, seed_cv) {
  n_full <- nrow(x_full)
  has_seed <- !is.null(seed_cv)
  if (has_seed) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(as.integer(seed_cv) + 1009L)   # same offset as bagging loop
  }
  for (b in seq_along(fits)) {
    idx <- sample.int(n_full, n_full, replace = TRUE)
    x_b <- x_full[idx, , drop = FALSE]
    y_b <- y_full[idx]
    bx_b <- mars_basis_cpp(x_b, fits[[b]]$dirs, fits[[b]]$cuts,
                           as.integer(fits[[b]]$selected.terms))
    g <- tryCatch(
      stats::glm.fit(x = bx_b, y = y_b,
                     family = stats::binomial(), intercept = FALSE),
      error = function(e) NULL
    )
    if (!is.null(g)) fits[[b]]$coefficients <- as.numeric(g$coefficients)
  }
  fits
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
                         cv_1se = FALSE) {
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

  # Run one fold and return per-size MSE on holdout.
  run_fold <- function(train_idx, test_idx) {
    x_tr <- x[train_idx, , drop = FALSE]
    y_tr <- y[train_idx]
    x_te <- x[test_idx,  , drop = FALSE]
    y_te <- y[test_idx]
    # Note: nprune passed through. force_size = 0 (we want the full path);
    # return_path = 1 to get path.subsets / path.coefs.
    fit <- mars_fit_cpp(x_tr, as.numeric(y_tr), degree,
                        as.integer(nk %||% (min(200L, max(20L, 2L * ncol(x_tr))) + 1L)),
                        as.numeric(penalty %||% (if (degree > 1L) 3 else 2)),
                        thresh, minspan, endspan, adjust_endspan, auto_linpreds,
                        fast_k, fast_beta, as.integer(nprune %||% 0L),
                        0L, trace, nthreads, 0L, 1L)
    # path.subsets is length M; entries with non-empty subset correspond to
    # achieved sizes. Build bx_test once for each size and accumulate MSE.
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
      out_mse[s] <- mean(r * r)
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
                           as.integer(size_star), 1L)
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
                               as.integer(size_star), 1L)
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
                           warmstart = TRUE) {
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
        warmstart = FALSE
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
      0L, trace, nthreads, 0L, 0L
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
  score_one <- function(d, pen, nk_c, fk, tr, te) {
    fit_k <- mars_fit_cpp(
      x[tr, , drop = FALSE], as.numeric(y[tr]),
      as.integer(d), as.integer(nk_c), as.numeric(pen),
      thresh, minspan, endspan, adjust_endspan, auto_linpreds,
      as.integer(fk), fast_beta,
      if (length(nprune) == 0L || is.na(nprune[1])) as.integer(nk_c) else as.integer(nprune),
      0L, trace, nthreads, 0L, 0L
    )
    bx_te <- mars_basis_cpp(x[te, , drop = FALSE],
                            fit_k$dirs, fit_k$cuts,
                            as.integer(fit_k$selected.terms))
    yhat <- drop(bx_te %*% as.numeric(fit_k$coefficients))
    r_te <- y[te] - yhat
    mean(r_te * r_te)
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
          trace, nthreads, 0L, 0L
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

      # Per-cell backward replay with its own penalty.
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
            nthreads, 0L, 0L
          ),
          error = function(e) NULL
        )
        if (is.null(bk)) { fold_mse[mi, fp_i] <- NA_real_; next }
        bx_te <- mars_basis_cpp(x_te, dirs_full, cuts_full,
                                as.integer(bk$selected.terms))
        yhat <- drop(bx_te %*% as.numeric(bk$coefficients))
        r_te <- y_te - yhat
        fold_mse[mi, fp_i] <- mean(r_te * r_te)
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
    0L, trace, nthreads, 0L, 0L
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

