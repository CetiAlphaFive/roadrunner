# ares — main fitting function

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
#'   hinge would deepen an existing interaction (parent term has degree
#'   >= 1). Default `1` (no adjustment). Earth's analogous default is `2`;
#'   ares's heuristic is empirically mixed on the inst/sims grid (helps a
#'   few cells, hurts others), so the default is conservative. Pass `2`
#'   to enable.
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
#' @param trace Trace level. 0 = silent (default), 1 = forward-pass progress.
#' @param nthreads Number of threads. 0 (default) selects
#'   `RcppParallel::defaultNumThreads()`. CRAN-distributed examples and the
#'   bundled vignette cap this at 2.
#' @param ... Additional arguments. Currently ignored.
#' @return An object of class `"ares"` — a list with components
#'   `coefficients`, `bx`, `dirs`, `cuts`, `selected.terms`, `rss`, `gcv`,
#'   `rss.per.subset`, `gcv.per.subset`, `fitted.values`, `residuals`, `namesx`,
#'   `call`, plus echoed control parameters.
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics* 19(1):1-67.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' print(fit)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
#' @export
ares <- function(x, ...) UseMethod("ares")

#' @rdname ares
#' @export
ares.formula <- function(x, data = NULL, ..., y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()
  mf <- stats::model.frame(formula, data = data)
  yv <- stats::model.response(mf)
  if (is.null(yv)) stop("ares: response variable is missing from formula/data.")
  mm <- stats::model.matrix(formula, mf)
  # drop intercept column
  has_int <- "(Intercept)" %in% colnames(mm)
  if (has_int) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  out <- ares.default(x = mm, y = as.numeric(yv), ...)
  out$call <- cl
  out$terms <- stats::terms(formula, data = data)
  out
}

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
                         trace = 0L, nthreads = 0L, ...) {
  cl <- match.call()
  pmethod <- match.arg(pmethod)
  autotune_speed <- match.arg(autotune.speed)
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
  if (is.data.frame(x)) {
    nm <- names(x)
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    colnames(x) <- nm
  }
  if (!is.matrix(x) || !is.numeric(x))
    stop("ares: x must be a numeric matrix or data frame.")
  if (any(!is.finite(x)))
    stop("ares: x contains NA, NaN, or non-finite values.")
  if (!is.numeric(y))
    stop("ares: y must be numeric.")
  if (any(!is.finite(y)))
    stop("ares: y contains NA, NaN, or non-finite values.")
  if (length(y) != nrow(x))
    stop("ares: length(y) (", length(y), ") must equal nrow(x) (", nrow(x), ").")
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

  nthreads <- as.integer(nthreads)
  nthreads_eff <- if (nthreads <= 0L) RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  # ---- Capture column names ----
  namesx <- colnames(x)
  if (is.null(namesx)) namesx <- paste0("V", seq_len(p))

  # ---- Coerce x to plain numeric matrix ----
  storage.mode(x) <- "double"

  # ---- Autotune dispatch (Phase 2 — v0.15+) ----
  if (isTRUE(autotune)) {
    # Promote nfold to a sensible default if user didn't ask.
    nfold_at <- if (nfold > 0L) nfold else 5L
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
      autotune_speed = autotune_speed
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
  class(out) <- c("ares")
  out
}

# Internal: build term labels matching earth's "h(x-cut)" format
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
#   - Stratification: regression-only — quantile-bin y into nfold bins and
#     round-robin within bins. ncross > 1 reseeds the partition.
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
    stop("ares: CV produced no evaluations — check nfold / ncross.")
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
.ares_autotune <- function(x, y, nk, thresh, minspan, endspan,
                           adjust_endspan, auto_linpreds, fast_k, fast_beta,
                           nprune, nfold, ncross, stratify, seed_cv,
                           cv_1se, trace, nthreads,
                           autotune_speed = "balanced") {
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
  nk_grid <- unique(pmin(c(nk_eff, 2L * nk_eff, 4L * nk_eff), 200L))

  # autotune_speed determines fast.k policy:
  #   "quality"  : fast.k = 0 always (no cache, slowest, most accurate).
  #   "fast"     : fast.k = 5 always (aggressive cache, cheapest).
  #   "balanced" : sweep fast.k in {10, 25, 0 (== "infinity" / no cache)}
  #                inside the grid; pick smallest within 1% of best.
  fk_grid <- switch(autotune_speed,
                    quality  = 0L,
                    fast     = 5L,
                    balanced = c(10L, 25L, 0L))

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

  # Successive-halving: score every cell on fold 1, build a running best;
  # any cell whose fold-1 MSE exceeds halving_factor * running_best is
  # eliminated and gets cv_mse = Inf (and no further folds). The halving
  # factor of 1.5 is empirical (loose enough not to drop close runners,
  # tight enough to drop wildly bad combos like deg=1 on a clear
  # interaction DGP).
  halving_factor <- 1.5
  ncell <- length(cells)
  fold_pairs <- list()  # ordered list of (rep_idx, k) folds
  for (r in seq_len(ncross)) {
    for (k in seq_len(nfold)) {
      fold_pairs[[length(fold_pairs) + 1L]] <- c(r, k)
    }
  }
  fold_mse <- matrix(NA_real_, nrow = ncell, ncol = length(fold_pairs))
  alive <- rep(TRUE, ncell)
  for (fp_i in seq_along(fold_pairs)) {
    r_idx <- fold_pairs[[fp_i]][1]
    k_idx <- fold_pairs[[fp_i]][2]
    fid <- fold_lists[[r_idx]]
    tr <- which(fid != k_idx); te <- which(fid == k_idx)
    if (length(te) < 1L) next
    for (i in seq_len(ncell)) {
      if (!alive[i]) next
      ce <- cells[[i]]
      m <- tryCatch(
        score_one(ce$degree, ce$penalty, ce$nk, ce$fast_k, tr, te),
        error = function(e) NA_real_
      )
      fold_mse[i, fp_i] <- m
    }
    # Halving check after fold 1 only (fp_i == 1).
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
                     vapply(cells, function(c) c$degree, integer(1L)),
                     vapply(cells, function(c) c$nk, integer(1L)),
                     vapply(cells, function(c) c$penalty, numeric(1L)),
                     vapply(cells, function(c) c$fast_k, integer(1L)))

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
      fk_vec <- vapply(cells, function(c) c$fast_k, integer(1L))
      score_key <- ifelse(fk_vec == 0L, .Machine$integer.max, fk_vec)
      candidates <- which(within1)
      ord_within <- candidates[order(score_key[candidates],
                                     scores[candidates],
                                     vapply(cells[candidates],
                                            function(c) c$degree,
                                            integer(1L)),
                                     vapply(cells[candidates],
                                            function(c) c$nk,
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
    degree     = vapply(cells, function(c) c$degree, integer(1L)),
    penalty    = vapply(cells, function(c) c$penalty, numeric(1L)),
    nk         = vapply(cells, function(c) c$nk, integer(1L)),
    fast_k     = vapply(cells, function(c) c$fast_k, integer(1L)),
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
    speed         = autotune_speed
  )
  fit_full
}

