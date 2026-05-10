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
                         seed.cv = NULL,
                         trace = 0L, nthreads = 0L, ...) {
  cl <- match.call()
  pmethod <- match.arg(pmethod)
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

  # ---- Call C++ engine ----
  if (pmethod == "cv") {
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
                        seed_cv = seed.cv)
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
  out$pmethod <- pmethod
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
                         nfold, ncross, stratify, seed_cv) {
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
  # Prefer sizes achieved by ALL folds (otherwise mean is biased toward
  # folds whose forward pass happened to keep more terms). Fall back to
  # any-fold sizes only if no size was achieved by every fold.
  full_n <- nrow(mse_mat)
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
  size_star <- which.min(cv_mse_finite)

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
    nfold      = nfold,
    ncross     = ncross,
    stratify   = stratify,
    cv.mse     = cv_mse,
    cv.n       = cv_n,
    size.star  = size_star
  )
  fit_full
}

# Tiny helper: %||% (used only inside .ares_cv_fit for safe defaulting).
`%||%` <- function(a, b) if (is.null(a)) b else a

