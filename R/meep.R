# ============================================================================
#  meep() -- cross-fitted causal ensemble (Phase Q-DML)
# ============================================================================
#
# Pure-R orchestration layer over the package's `ares()` (MARS) and `krls()`
# (Kernel Regularized Least Squares) fitters. `meep()` produces cross-fitted
# (out-of-fold) nuisance estimates suitable for Double Machine Learning and
# causal forests -- it does NOT estimate treatment effects itself. The honest
# OOF prediction matrix it returns is meant to be handed downstream to
# `grf::causal_forest()` (via `Y.hat` / `W.hat`) or to a DML routine.
#
# No C++ in this file. No new hard dependencies: the non-negative least
# squares solver used by the stacking ensemble is hand-rolled below
# (Lawson-Hanson active-set).
#
# File contents:
#   - meep()                 main entry point
#   - build_folds()          K-fold, cluster-aware, seeded fold construction
#   - .meep_nnls()           Lawson-Hanson non-negative least squares
#   - .meep_learner_specs()  named-list dispatch of adapter closures
#   - .meep_combine()        stacking / averaging / best ensemble combine
#   - predict.meep()         predictions on NEW data (full-data refits)
#   - print.meep()           compact summary
#   - summary.meep()         full diagnostics incl. overlap
#   - print.summary.meep()   formatter

# ----------------------------------------------------------------------------
#  Small internal utilities
# ----------------------------------------------------------------------------
# (`%||%` is defined in R/ares.R and shared package-wide.)

# Is a vector binary (exactly two distinct finite values, coercible to 0/1)?
.meep_is_binary <- function(v) {
  u <- unique(v[is.finite(v)])
  length(u) == 2L
}

# Coerce a 2-valued vector to 0/1 numeric. The smaller value maps to 0.
.meep_to01 <- function(v) {
  u <- sort(unique(v[is.finite(v)]))
  as.numeric(v == u[2L])
}

# Is learner `name` applicable to a nuisance of family `fam`?
#
# `ares` and `krls` are family-agnostic (always apply). `ols` is a gaussian
# learner only; `logreg` is a binomial learner only. Any other name is a
# custom list-learner -- not filtered (the user is responsible for matching
# it to the nuisance families).
.meep_learner_family_ok <- function(name, fam) {
  if (identical(name, "ols"))    return(identical(fam, "gaussian"))
  if (identical(name, "logreg")) return(identical(fam, "binomial"))
  # ares, krls, and any custom list-learner: applicable to both families.
  TRUE
}

# ----------------------------------------------------------------------------
#  build_folds() -- K-fold assignment, cluster-aware, seeded
# ----------------------------------------------------------------------------

#' Construct cross-fitting fold assignments
#'
#' Builds a length-`n` integer vector of fold ids in `1:K`. When `cluster`
#' is supplied, whole clusters are kept together within a fold (folds are
#' assigned at the cluster level), which is required for cluster-robust
#' cross-fitting. Used internally by [meep()]; exported is unnecessary.
#'
#' @param n Integer sample size.
#' @param K Integer number of folds.
#' @param cluster Optional length-`n` grouping vector; clusters kept whole.
#' @param seed Optional integer seed for reproducible fold draws.
#' @return Integer vector of length `n` with values in `1:K`.
#' @keywords internal
#' @noRd
build_folds <- function(n, K, cluster = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- as.integer(K)
  if (is.na(K) || K < 2L)
    stop("meep: `folds` must be an integer >= 2 (got ", K, ").",
         call. = FALSE)
  if (K > n)
    stop("meep: `folds` (", K, ") cannot exceed the sample size (", n, ").",
         call. = FALSE)

  if (is.null(cluster)) {
    # plain K-fold: shuffle then deal round-robin so fold sizes are balanced
    ord <- sample.int(n)
    fid <- integer(n)
    fid[ord] <- rep_len(seq_len(K), n)
    return(fid)
  }

  if (length(cluster) != n)
    stop("meep: `cluster` must have length n = ", n, " (got ",
         length(cluster), ").", call. = FALSE)
  cl <- as.character(cluster)
  ucl <- unique(cl)
  ncl <- length(ucl)
  if (ncl < K)
    stop("meep: number of clusters (", ncl,
         ") is smaller than `folds` (", K, ").", call. = FALSE)
  # assign each cluster to a fold, round-robin over a shuffled cluster order
  cl_ord <- sample.int(ncl)
  cl_fold <- integer(ncl)
  cl_fold[cl_ord] <- rep_len(seq_len(K), ncl)
  names(cl_fold) <- ucl
  as.integer(cl_fold[cl])
}

# ----------------------------------------------------------------------------
#  .meep_nnls() -- non-negative least squares (Lawson-Hanson active set)
# ----------------------------------------------------------------------------
#
# Solves   min_w ||A w - b||^2   s.t.   w >= 0
# Hand-rolled, no dependencies. Standard Lawson & Hanson (1974) algorithm.

.meep_nnls <- function(A, b, tol = 1e-10, maxit = NULL) {
  A <- as.matrix(A)
  m <- nrow(A); p <- ncol(A)
  if (is.null(maxit)) maxit <- 3L * p
  w <- numeric(p)                       # solution
  P <- logical(p)                       # passive (free) set membership
  AtA <- crossprod(A)
  Atb <- as.numeric(crossprod(A, b))
  grad <- Atb - as.numeric(AtA %*% w)   # negative gradient of 0.5||Aw-b||^2
  iter <- 0L
  while (any(!P) && max(grad[!P]) > tol && iter < 30L * maxit) {
    iter <- iter + 1L
    # move the most-violating variable into the passive set
    cand <- which(!P)
    j <- cand[which.max(grad[cand])]
    P[j] <- TRUE
    # inner loop: solve unconstrained LS on the passive set
    repeat {
      idx <- which(P)
      z <- numeric(p)
      sol <- tryCatch(
        solve(AtA[idx, idx, drop = FALSE], Atb[idx]),
        error = function(e) NULL)
      if (is.null(sol)) {
        # singular passive subsystem: fall back to a tiny ridge
        rg <- diag(1e-8, length(idx))
        sol <- solve(AtA[idx, idx, drop = FALSE] + rg, Atb[idx])
      }
      z[idx] <- sol
      if (all(z[idx] > tol)) {
        w <- z
        break
      }
      # some passive coefficient went non-positive: step back
      neg <- idx[z[idx] <= tol]
      alpha <- min(w[neg] / (w[neg] - z[neg]))
      w <- w + alpha * (z - w)
      P[abs(w) < tol] <- FALSE
      if (!any(P)) { w <- z; break }
    }
    grad <- Atb - as.numeric(AtA %*% w)
  }
  w[w < 0] <- 0
  w
}

# ----------------------------------------------------------------------------
#  Learner adapter dispatch
# ----------------------------------------------------------------------------
#
# Each learner is a named entry with a $fit closure and a $predict closure.
# The named-list design is the extensibility hook: a future "ranger" entry
# slots in here with no change to the cross-fitting machinery.
#
#   $fit(X, y, family, weights, hp)  -> a fitted model object
#   $predict(model, newX, family)    -> numeric vector (response scale;
#                                       probabilities for binomial)

#
# `ols` and `logreg` are opt-in (not in the default `learners`): they are
# unregularized linear fitters with NO hyperparameters, so `tune` is a
# no-op for them -- every fold is a plain refit. `ols` is the gaussian
# learner; `logreg` is the binomial learner (binary outcome / propensity).
.meep_learner_specs <- function(ares_args = list(), krls_args = list(),
                                ols_args = list(), logreg_args = list()) {
  list(
    ares = list(
      fit = function(X, y, family, weights, hp) {
        args <- c(list(x = X, y = y, family = family),
                  if (!is.null(weights)) list(weights = weights),
                  ares_args, hp)
        do.call(ares, args)
      },
      predict = function(model, newX, family) {
        as.numeric(predict(model, newdata = newX, type = "response"))
      }
    ),
    krls = list(
      fit = function(X, y, family, weights, hp) {
        loss <- if (identical(family, "binomial")) "logistic" else "ls"
        args <- c(list(X = X, y = y, loss = loss),
                  if (!is.null(weights)) list(weights = weights),
                  krls_args, hp)
        # krls is verbose-by-default on some paths; keep it quiet here
        do.call(krls, args)
      },
      predict = function(model, newX, family) {
        ty <- if (identical(family, "binomial")) "prob" else "response"
        pr <- predict(model, newdata = newX, type = ty)
        # predict.krls_rr returns a list with a $fit matrix; predict.ares
        # returns a bare numeric. Normalise to a numeric vector.
        if (is.list(pr)) pr <- pr$fit
        as.numeric(pr)
      }
    ),
    ols = list(
      # Gaussian learner. No hyperparameters: `hp` is always empty.
      fit = function(X, y, family, weights, hp) {
        if (identical(family, "binomial"))
          stop("meep: learner 'ols' cannot fit a binomial nuisance; ",
               "use 'logreg' for binary outcomes/treatments.",
               call. = FALSE)
        args <- c(list(x = X, y = y),
                  if (!is.null(weights)) list(weights = weights),
                  ols_args)
        do.call(ols, args)
      },
      predict = function(model, newX, family) {
        as.numeric(predict(model, newdata = newX, type = "response"))
      }
    ),
    logreg = list(
      # Binomial learner. No hyperparameters: `hp` is always empty.
      fit = function(X, y, family, weights, hp) {
        if (!identical(family, "binomial"))
          stop("meep: learner 'logreg' only fits binomial nuisances; ",
               "use 'ols' for a continuous outcome/treatment.",
               call. = FALSE)
        args <- c(list(x = X, y = y),
                  if (!is.null(weights)) list(weights = weights),
                  logreg_args)
        do.call(logreg, args)
      },
      predict = function(model, newX, family) {
        as.numeric(predict(model, newdata = newX, type = "response"))
      }
    )
  )
}

# ----------------------------------------------------------------------------
#  .meep_combine() -- ensemble combination over the OOF matrix
# ----------------------------------------------------------------------------
#
# OOF is an n x L matrix; column l is the honest cross-fitted prediction of
# learner l. Some columns may be entirely NA for some rows (a learner that
# failed a fold). Returns a named numeric weight vector summing to 1.

.meep_combine <- function(OOF, resp, ensemble, weights, cv_perf_mse) {
  L <- ncol(OOF)
  lnames <- colnames(OOF)
  # a learner is "alive" if it produced any non-NA prediction at all
  alive <- vapply(seq_len(L), function(l) any(is.finite(OOF[, l])),
                  logical(1L))
  if (!any(alive))
    stop("meep: every learner failed on every fold; cannot ensemble.",
         call. = FALSE)

  if (identical(ensemble, "average")) {
    w <- as.numeric(alive)
    w <- w / sum(w)
    names(w) <- lnames
    return(w)
  }

  if (identical(ensemble, "best")) {
    w <- numeric(L)
    perf <- cv_perf_mse
    perf[!alive] <- Inf
    w[which.min(perf)] <- 1
    names(w) <- lnames
    return(w)
  }

  # ensemble == "stack": NNLS of resp on the OOF columns, on rows where all
  # alive learners have a prediction (complete cases over the alive set).
  alive_idx <- which(alive)
  A <- OOF[, alive_idx, drop = FALSE]
  cc <- stats::complete.cases(A) & is.finite(resp)
  if (sum(cc) < 2L) {
    # not enough complete rows to stack: fall back to equal weights
    w <- as.numeric(alive)
    w <- w / sum(w)
    names(w) <- lnames
    return(w)
  }
  Acc <- A[cc, , drop = FALSE]
  bcc <- resp[cc]
  if (!is.null(weights)) {
    sw <- sqrt(weights[cc])
    Acc <- Acc * sw
    bcc <- bcc * sw
  }
  wa <- .meep_nnls(Acc, bcc)
  s <- sum(wa)
  if (s <= 0) {
    # degenerate: NNLS zeroed everything -> equal weight over alive learners
    wa <- rep(1, length(alive_idx))
    s <- sum(wa)
  }
  wa <- wa / s
  w <- numeric(L)
  w[alive_idx] <- wa
  names(w) <- lnames
  w
}

# Row-wise ensemble prediction: OOF %*% w, but renormalizing the weights
# over the non-NA columns for each row independently (graceful degradation).
.meep_apply_weights <- function(OOF, w) {
  n <- nrow(OOF); L <- ncol(OOF)
  out <- numeric(n)
  for (i in seq_len(n)) {
    row <- OOF[i, ]
    ok <- is.finite(row) & w > 0
    if (!any(ok)) {
      # no learner predicted this row at all -- should not happen because
      # an all-fail fold stops earlier, but guard anyway
      out[i] <- NA_real_
      next
    }
    wi <- w[ok] / sum(w[ok])
    out[i] <- sum(row[ok] * wi)
  }
  out
}

# ----------------------------------------------------------------------------
#  meep() -- main entry point
# ----------------------------------------------------------------------------

#' Cross-fitted causal ensemble of `ares()` and `krls()`
#'
#' `meep()` produces honest, cross-fitted (out-of-fold) nuisance estimates
#' for use in Double Machine Learning (DML) and causal-forest workflows.
#' For an outcome `y`, an optional `treatment`, and covariates `X`, it
#' cross-fits an ensemble of the package's [ares()] (MARS) and [krls()]
#' (Kernel Regularized Least Squares) learners and returns the
#' out-of-fold predictions \eqn{\hat E[Y\mid X]} and (when `treatment` is
#' supplied) \eqn{\hat E[D\mid X]}, plus the residuals \eqn{Y-\hat Y} and
#' \eqn{D-\hat D}.
#'
#' `meep()` does **not** estimate a treatment effect. The honest OOF
#' predictions are meant to be handed to a downstream DML estimator or to
#' `grf::causal_forest()` (see Examples). Because each column of the OOF
#' prediction matrix already comes from a model that never saw the row it
#' predicts, the stacking ensemble (`ensemble = "stack"`) can fit
#' non-negative least squares directly on that matrix without a second
#' layer of nesting -- leakage is already purged column-wise.
#'
#' @section Family auto-detection:
#' When `family` is `NULL`, a binary `y` (exactly two distinct values)
#' triggers `ares(family = "binomial")` and `krls(loss = "logistic")`;
#' otherwise the gaussian path is used. The same logic applies to
#' `treatment_family`: a binary treatment is modelled as a propensity
#' score via logistic regression, a continuous treatment via gaussian
#' regression. A treatment with three or more distinct values is rejected
#' (multi-valued / continuous-treatment IRM is out of scope).
#'
#' The family also governs which default learners are fitted for each
#' nuisance. `ares` and `krls` are family-agnostic and always apply.
#' `ols` is fitted only for gaussian (continuous) nuisances and `logreg`
#' only for binomial (binary) nuisances. Because the Y-model and the
#' D-model can differ in family, applicability is resolved *per nuisance*:
#' a gaussian outcome with a binary treatment fits `ols` for
#' `outcome`/`mu0`/`mu1` and `logreg` for `treatment`. A learner that
#' does not apply to a nuisance is simply skipped -- its OOF column stays
#' all-`NA` and it is *not* recorded as a fold failure.
#'
#' @section Tuning (`tune`):
#' * `"once"` (default) -- autotune each learner on the full data, freeze
#'   the resulting hyperparameters, and refit those frozen settings within
#'   every fold. Standard DML practice: fit honesty (cross-fitting) is what
#'   protects inference; freezing hyperparameters is a deliberate, second-
#'   order compromise.
#' * `"per_fold"` -- autotune independently inside every fold. Purist
#'   option; slower.
#' * `"none"` -- the fast path. Call [ares()] and [krls()] with their own
#'   default arguments (plus anything passed via `ares_args` / `krls_args`).
#'   No autotune is invoked at all.
#'
#' @section Graceful degradation:
#' If a learner errors on a fold (for example, a constant nuisance in a
#' fold subset), its column for that fold becomes `NA`, the ensemble
#' weights are renormalized over the surviving learners on a per-row basis,
#' and the event is logged in `$fold_failures`. If *every* learner fails a
#' fold, `meep()` stops.
#'
#' @section DoubleML recipe:
#' To use these nuisances with the \pkg{DoubleML} package, pass
#' `y_hat_oof` / `d_hat_oof` as `external_predictions` and reuse
#' `fit$folds` as the sample split so the cross-fitting partition matches.
#'
#' @param X An n-by-p numeric matrix or data.frame of covariates. Factor
#'   columns in a data.frame are expanded by the base fitters.
#' @param y Length-n outcome vector.
#' @param treatment Optional length-n treatment / exposure vector. `NULL`
#'   gives an outcome-only fit (no `d_hat_oof`, no arm models).
#' @param folds Either a single integer `K` (number of folds), or a
#'   length-n integer vector of fold ids that is honored verbatim.
#' @param learners Character subset of `c("ares", "krls", "ols",
#'   "logreg")`, or a named list of adapter closures (extensibility
#'   hook). The default is all four: `c("ares", "krls", "ols",
#'   "logreg")`. `ares` and `krls` are family-agnostic and apply to
#'   every nuisance. `"ols"` and `"logreg"` are unregularized linear
#'   learners with no hyperparameters (so `tune` is a no-op for them --
#'   a plain refit per fold); they are applied *per nuisance by family*,
#'   `ols` for gaussian (continuous) targets and `logreg` for binomial
#'   (binary) ones. A learner that does not apply to a nuisance is
#'   skipped, leaving an all-`NA` OOF column that the ensemble excludes;
#'   it is not counted as a fold failure. Narrowing `learners` to names
#'   that are all family-incompatible with a nuisance (for example
#'   `learners = "ols"` with a binary outcome) is an error. Custom
#'   list-learners are never family-filtered.
#' @param ensemble How to combine learners: `"stack"` (non-negative least
#'   squares on the OOF matrix), `"average"` (equal weights), or `"best"`
#'   (pick the single lowest-OOF-loss learner).
#' @param arm_models `"auto"` fits arm-specific outcome models
#'   (\eqn{\mu_0}, \eqn{\mu_1}) only when the treatment is binary;
#'   `"always"` forces them; `"never"` suppresses them.
#' @param family Outcome family, `"gaussian"` or `"binomial"`. `NULL`
#'   auto-detects.
#' @param treatment_family Treatment family, `"gaussian"` or `"binomial"`.
#'   `NULL` auto-detects.
#' @param tune One of `"once"`, `"per_fold"`, `"none"`. See *Tuning*.
#' @param cluster Optional length-n grouping vector. Folds are assigned at
#'   the cluster level so each cluster lands entirely within one fold.
#' @param weights Optional length-n case weights, forwarded to every base
#'   fit and into the NNLS stacking objective.
#' @param seed Optional integer seed for fold construction.
#' @param ares_args A named list of extra arguments spliced into every
#'   [ares()] call.
#' @param krls_args A named list of extra arguments spliced into every
#'   [krls()] call.
#' @param ols_args A named list of extra arguments spliced into every
#'   [ols()] call (only relevant when `"ols"` is among `learners`).
#' @param logreg_args A named list of extra arguments spliced into every
#'   [logreg()] call (only relevant when `"logreg"` is among
#'   `learners`).
#' @param verbose Logical; if `TRUE`, print progress per nuisance / fold.
#' @param ... Currently unused; reserved.
#'
#' @return An object of S3 class `"meep"`: a list with `y_hat_oof`,
#'   `d_hat_oof`, `mu0_hat_oof`, `mu1_hat_oof`, `y_resid`, `d_resid`,
#'   `folds`, `oof_matrix` (per-nuisance n-by-L matrices), `ensemble_weights`,
#'   `learner_cv_perf`, `learners`, `family`, `treatment_family`, `tune`,
#'   `ensemble`, `n`, `p`, `frozen_hyperparams`, `fold_failures`,
#'   `cluster`, `weights`, `seed`, and `call`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 800; p <- 5
#' X <- matrix(runif(n * p), n, p)
#' m  <- 0.4 * X[, 1] + 0.3 * sin(3 * X[, 2])
#' D  <- m + rnorm(n, sd = 0.3)
#' g  <- X[, 1]^2 + 0.5 * X[, 3]
#' Y  <- 1.0 * D + g + rnorm(n, sd = 0.5)        # true ATE = 1.0
#' fit <- meep(X, Y, treatment = D, folds = 5, seed = 42)
#' fit
#' # Hand the honest OOF nuisances to a causal forest:
#' # grf::causal_forest(X, Y, D, Y.hat = fit$y_hat_oof, W.hat = fit$d_hat_oof)
#' }
#'
#' @seealso [ares()], [krls()]
#' @export
meep <- function(X, y, treatment = NULL,
                 folds = 5L,
                 learners = c("ares", "krls", "ols", "logreg"),
                 ensemble = c("stack", "average", "best"),
                 arm_models = c("auto", "always", "never"),
                 family = NULL,
                 treatment_family = NULL,
                 tune = c("once", "per_fold", "none"),
                 cluster = NULL,
                 weights = NULL,
                 seed = NULL,
                 ares_args = list(),
                 krls_args = list(),
                 ols_args = list(),
                 logreg_args = list(),
                 verbose = FALSE, ...) {
  cl <- match.call()
  ensemble   <- match.arg(ensemble)
  arm_models <- match.arg(arm_models)
  tune       <- match.arg(tune)

  # ---- input validation --------------------------------------------------
  if (is.data.frame(X)) {
    Xn <- nrow(X)
  } else {
    X <- as.matrix(X)
    if (!is.numeric(X))
      stop("meep: `X` must be numeric (or a data.frame).", call. = FALSE)
    Xn <- nrow(X)
  }
  n <- Xn
  p <- ncol(X)
  if (length(y) != n)
    stop("meep: length(y) (", length(y), ") must equal nrow(X) (", n, ").",
         call. = FALSE)
  y <- as.numeric(y)

  if (!is.null(treatment)) {
    if (length(treatment) != n)
      stop("meep: length(treatment) (", length(treatment),
           ") must equal nrow(X) (", n, ").", call. = FALSE)
  }
  if (!is.null(weights)) {
    if (length(weights) != n)
      stop("meep: length(weights) must equal nrow(X).", call. = FALSE)
    weights <- as.numeric(weights)
  }

  # ---- resolve learner specs --------------------------------------------
  if (is.list(learners)) {
    learner_specs <- learners
    learner_names <- names(learner_specs)
    if (is.null(learner_names) || any(!nzchar(learner_names)))
      stop("meep: a list of `learners` must be fully named.", call. = FALSE)
  } else {
    learners <- match.arg(learners, c("ares", "krls", "ols", "logreg"),
                          several.ok = TRUE)
    learners <- unique(learners)
    all_specs <- .meep_learner_specs(ares_args, krls_args,
                                     ols_args, logreg_args)
    learner_specs <- all_specs[learners]
    learner_names <- learners
  }
  L <- length(learner_specs)
  if (L < 1L)
    stop("meep: at least one learner is required.", call. = FALSE)

  # ---- resolve outcome family -------------------------------------------
  if (is.null(family)) {
    family <- if (.meep_is_binary(y)) "binomial" else "gaussian"
  } else {
    family <- match.arg(family, c("gaussian", "binomial"))
  }
  y_model <- if (identical(family, "binomial")) .meep_to01(y) else y

  # ---- resolve treatment family -----------------------------------------
  d_binary <- FALSE
  if (!is.null(treatment)) {
    u_treat <- unique(treatment[is.finite(treatment)])
    if (length(u_treat) < 2L)
      stop("meep: `treatment` is constant; nothing to model.", call. = FALSE)
    if (is.null(treatment_family)) {
      if (length(u_treat) == 2L) {
        treatment_family <- "binomial"
      } else if (is.factor(treatment) || is.character(treatment)) {
        stop("meep: `treatment` has ", length(u_treat),
             " distinct values. Multi-valued / categorical treatments are ",
             "not supported (only binary or continuous). ",
             "See the package scope notes.", call. = FALSE)
      } else if (length(u_treat) <= 10L &&
                 all(u_treat == round(u_treat))) {
        # a small set of integer-valued treatment levels is treated as
        # multi-valued and rejected (continuous treatment should have many
        # distinct values).
        stop("meep: `treatment` has ", length(u_treat),
             " distinct integer values. Multi-valued treatments are not ",
             "supported (only binary or continuous). For a genuinely ",
             "continuous treatment pass it as a numeric vector with many ",
             "distinct values.", call. = FALSE)
      } else {
        treatment_family <- "gaussian"
      }
    } else {
      treatment_family <- match.arg(treatment_family,
                                    c("gaussian", "binomial"))
    }
    d_binary <- identical(treatment_family, "binomial")
    treatment <- as.numeric(treatment)
    d_model <- if (d_binary) .meep_to01(treatment) else treatment
  } else {
    d_model <- NULL
  }

  # ---- decide whether arm models are engaged ----------------------------
  use_arms <- switch(arm_models,
                     auto   = !is.null(treatment) && d_binary,
                     always = !is.null(treatment),
                     never  = FALSE)
  if (identical(arm_models, "always") && !is.null(treatment) && !d_binary)
    stop("meep: arm_models = 'always' requires a binary treatment.",
         call. = FALSE)

  # ---- fold construction -------------------------------------------------
  if (length(folds) == n) {
    fold_id <- as.integer(folds)
    if (anyNA(fold_id))
      stop("meep: user-supplied `folds` vector contains NA.", call. = FALSE)
    K <- length(unique(fold_id))
    if (K < 2L)
      stop("meep: user-supplied `folds` must contain >= 2 distinct folds.",
           call. = FALSE)
  } else {
    if (length(folds) != 1L)
      stop("meep: `folds` must be a single integer K or a length-n vector.",
           call. = FALSE)
    K <- as.integer(folds)
    fold_id <- build_folds(n, K, cluster = cluster, seed = seed)
  }
  fold_levels <- sort(unique(fold_id))
  K <- length(fold_levels)

  # tiny-fold sanity warning
  for (k in fold_levels) {
    ntr <- sum(fold_id != k)
    if (ntr < 2L * p)
      warning("meep: fold ", k, " leaves only ", ntr,
              " training rows for p = ", p,
              " predictors (< 2p); fits may be unstable.", call. = FALSE)
  }

  # ---- tune == "once": autotune on full data, freeze hyperparameters ----
  frozen <- NULL
  if (identical(tune, "once")) {
    frozen <- .meep_freeze_hyperparams(X, y_model, family, weights,
                                       learner_names, ares_args, krls_args,
                                       verbose)
  }

  # ---- assemble the nuisance task list ----------------------------------
  nuisances <- list(
    outcome = list(resp = y_model, fam = family, rows = NULL)
  )
  if (!is.null(treatment)) {
    nuisances$treatment <- list(resp = d_model, fam = treatment_family,
                                rows = NULL)
    if (use_arms) {
      # arm 0 outcome model: trained on rows with D == 0
      nuisances$mu0 <- list(resp = y_model, fam = family,
                            rows = which(d_model == 0))
      nuisances$mu1 <- list(resp = y_model, fam = family,
                            rows = which(d_model == 1))
    }
  }

  # ---- cross-fitting loop -----------------------------------------------
  oof_matrix       <- list()
  ensemble_weights <- list()
  oof_pred         <- list()
  fold_failures    <- list()
  cv_perf_rows     <- list()

  for (nm in names(nuisances)) {
    N <- nuisances[[nm]]
    resp <- N$resp
    fam  <- N$fam
    OOF  <- matrix(NA_real_, n, L,
                   dimnames = list(NULL, learner_names))
    # which rows are eligible to be predicted (arm models: only that arm's
    # rows have a defined nuisance; we still predict the whole fold so the
    # OOF vector is length n, but residuals/usage downstream is the user's
    # job -- arm models cross-fit on arm-restricted training rows).
    arm_rows <- N$rows

    # which learners apply to this nuisance's family. `ols` is gaussian-only,
    # `logreg` binomial-only; an inapplicable learner is skipped (its OOF
    # column stays all-NA) -- not a failure. Custom list-learners always
    # apply. An empty applicable set means the user narrowed `learners` to
    # family-incompatible names; stop with a clear message.
    applicable <- vapply(learner_names,
                         function(ln) .meep_learner_family_ok(ln, fam),
                         logical(1L))
    if (!any(applicable))
      stop("meep: no applicable learner for nuisance '", nm,
           "' (family '", fam, "'). The requested learners (",
           paste(learner_names, collapse = ", "),
           ") are all family-incompatible: 'ols' needs a continuous ",
           "(gaussian) target, 'logreg' a binary (binomial) one.",
           call. = FALSE)

    for (k in fold_levels) {
      te <- which(fold_id == k)
      tr <- which(fold_id != k)
      if (!is.null(arm_rows)) {
        tr <- intersect(tr, arm_rows)
      }
      if (length(tr) < 2L) {
        # not enough training rows for this fold/arm: whole fold is NA.
        # Only applicable learners are recorded as failures -- a learner
        # skipped for family incompatibility is not a failure.
        for (li in seq_len(L)) {
          if (!applicable[li]) next
          fold_failures[[length(fold_failures) + 1L]] <- list(
            nuisance = nm, fold = k, learner = learner_names[li],
            reason = "insufficient training rows")
        }
        next
      }
      Xtr <- if (is.data.frame(X)) X[tr, , drop = FALSE] else
        X[tr, , drop = FALSE]
      Xte <- if (is.data.frame(X)) X[te, , drop = FALSE] else
        X[te, , drop = FALSE]
      ytr <- resp[tr]
      wtr <- if (is.null(weights)) NULL else weights[tr]

      n_failed_here <- 0L
      for (li in seq_len(L)) {
        # skip a learner that does not apply to this nuisance's family;
        # its OOF column stays all-NA and no failure is recorded.
        if (!applicable[li]) next
        lname <- learner_names[li]
        spec  <- learner_specs[[li]]
        hp <- .meep_pick_hp(tune, frozen, lname,
                            Xtr, ytr, fam, wtr,
                            ares_args, krls_args)
        fit_l <- tryCatch(
          spec$fit(Xtr, ytr, fam, wtr, hp),
          error = function(e) e)
        if (inherits(fit_l, "error")) {
          fold_failures[[length(fold_failures) + 1L]] <- list(
            nuisance = nm, fold = k, learner = lname,
            reason = conditionMessage(fit_l))
          n_failed_here <- n_failed_here + 1L
          next
        }
        pr <- tryCatch(
          spec$predict(fit_l, Xte, fam),
          error = function(e) e)
        if (inherits(pr, "error")) {
          fold_failures[[length(fold_failures) + 1L]] <- list(
            nuisance = nm, fold = k, learner = lname,
            reason = paste0("predict: ", conditionMessage(pr)))
          n_failed_here <- n_failed_here + 1L
          next
        }
        OOF[te, li] <- as.numeric(pr)
        if (verbose)
          message("  [", nm, "] fold ", k, " / learner ", lname, " ok")
      }
      if (n_failed_here == sum(applicable))
        stop("meep: every learner failed on fold ", k,
             " for nuisance '", nm, "'. Cannot cross-fit.", call. = FALSE)
    }

    # per-learner OOF loss (MSE for gaussian, binomial deviance for binary)
    perf <- numeric(L)
    for (li in seq_len(L)) {
      col <- OOF[, li]
      ok <- is.finite(col) & is.finite(resp)
      if (!is.null(arm_rows)) ok <- ok & (seq_len(n) %in% arm_rows)
      if (!any(ok)) { perf[li] <- Inf; next }
      perf[li] <- .meep_loss(resp[ok], col[ok], fam,
                             if (is.null(weights)) NULL else weights[ok])
      cv_perf_rows[[length(cv_perf_rows) + 1L]] <- data.frame(
        nuisance = nm, learner = learner_names[li],
        loss = perf[li],
        metric = if (identical(fam, "binomial")) "deviance" else "mse",
        stringsAsFactors = FALSE)
    }

    # combine -- weights are scored on the rows relevant to this nuisance
    if (!is.null(arm_rows)) {
      OOF_score <- OOF[arm_rows, , drop = FALSE]
      resp_score <- resp[arm_rows]
      w_score <- if (is.null(weights)) NULL else weights[arm_rows]
    } else {
      OOF_score <- OOF
      resp_score <- resp
      w_score <- weights
    }
    wts <- .meep_combine(OOF_score, resp_score, ensemble, w_score, perf)
    pred <- .meep_apply_weights(OOF, wts)

    oof_matrix[[nm]]       <- OOF
    ensemble_weights[[nm]] <- wts
    oof_pred[[nm]]         <- pred
  }

  learner_cv_perf <- if (length(cv_perf_rows))
    do.call(rbind, cv_perf_rows) else
    data.frame(nuisance = character(0), learner = character(0),
               loss = numeric(0), metric = character(0),
               stringsAsFactors = FALSE)
  rownames(learner_cv_perf) <- NULL

  # ---- assemble return ---------------------------------------------------
  y_hat_oof <- oof_pred$outcome
  d_hat_oof <- oof_pred$treatment %||% NULL
  mu0_hat_oof <- oof_pred$mu0 %||% NULL
  mu1_hat_oof <- oof_pred$mu1 %||% NULL

  y_resid <- y_model - y_hat_oof
  d_resid <- if (!is.null(d_hat_oof)) d_model - d_hat_oof else NULL

  out <- list(
    y_hat_oof        = y_hat_oof,
    d_hat_oof        = d_hat_oof,
    mu0_hat_oof      = mu0_hat_oof,
    mu1_hat_oof      = mu1_hat_oof,
    y_resid          = y_resid,
    d_resid          = d_resid,
    folds            = fold_id,
    oof_matrix       = oof_matrix,
    ensemble_weights = ensemble_weights,
    learner_cv_perf  = learner_cv_perf,
    learners         = learner_names,
    family           = family,
    treatment_family = treatment_family,
    tune             = tune,
    ensemble         = ensemble,
    arm_models       = arm_models,
    n                = n,
    p                = p,
    frozen_hyperparams = frozen,
    fold_failures    = fold_failures,
    cluster          = cluster,
    weights          = weights,
    seed             = seed,
    call             = cl
  )
  # lazy cache env for predict.meep() full-data refits (not serialized)
  out$.full_fits <- new.env(parent = emptyenv())
  # stash data needed for full-data refits
  out$.train <- list(X = X, y = y_model, d = d_model,
                      ares_args = ares_args, krls_args = krls_args,
                      learner_specs = learner_specs)
  class(out) <- "meep"
  out
}

# ----------------------------------------------------------------------------
#  Hyperparameter freezing / picking helpers
# ----------------------------------------------------------------------------

# Autotune each learner once on the full data and capture the winning
# hyperparameters as a per-learner list of argument lists.
.meep_freeze_hyperparams <- function(X, y, family, weights, learner_names,
                                     ares_args, krls_args, verbose) {
  frozen <- list()
  for (lname in learner_names) {
    if (identical(lname, "ares")) {
      args <- c(list(x = X, y = y, family = family, autotune = TRUE),
                if (!is.null(weights)) list(weights = weights),
                ares_args)
      fit <- tryCatch(do.call(ares, args), error = function(e) e)
      if (inherits(fit, "error") || is.null(fit$autotune)) {
        frozen[[lname]] <- list()
      } else {
        at <- fit$autotune
        frozen[[lname]] <- list(degree = at$degree, penalty = at$penalty,
                                nk = at$nk, fast.k = at$fast_k)
      }
    } else if (identical(lname, "krls")) {
      loss <- if (identical(family, "binomial")) "logistic" else "ls"
      args <- c(list(X = X, y = y, loss = loss, autotune = TRUE),
                if (!is.null(weights)) list(weights = weights),
                krls_args)
      fit <- tryCatch(do.call(krls, args), error = function(e) e)
      if (inherits(fit, "error")) {
        frozen[[lname]] <- list()
      } else {
        sig <- fit$sigma_vec %||% fit$sigma
        frozen[[lname]] <- list(sigma = sig, lambda = fit$lambda)
      }
    } else {
      # ols / logreg (no hyperparameters) and any unknown learner --
      # nothing to freeze; the per-fold call is a plain refit.
      frozen[[lname]] <- list()
    }
    if (verbose)
      message("  [tune=once] froze hyperparameters for learner ", lname)
  }
  frozen
}

# Decide the hyperparameter list to splice into a per-fold base call.
.meep_pick_hp <- function(tune, frozen, lname, Xtr, ytr, fam, wtr,
                          ares_args, krls_args) {
  if (identical(tune, "none")) {
    # fast path: pure base-learner defaults. No autotune, no frozen hp.
    return(list())
  }
  if (identical(tune, "once")) {
    hp <- frozen[[lname]]
    if (is.null(hp)) hp <- list()
    return(hp)
  }
  # tune == "per_fold": autotune inside this fold
  if (identical(lname, "ares")) {
    return(list(autotune = TRUE))
  }
  if (identical(lname, "krls")) {
    return(list(autotune = TRUE))
  }
  list()
}

# ----------------------------------------------------------------------------
#  Loss helper
# ----------------------------------------------------------------------------

.meep_loss <- function(resp, pred, family, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(resp))
  if (identical(family, "binomial")) {
    eps <- 1e-12
    p <- pmin(pmax(pred, eps), 1 - eps)
    dev <- -2 * (resp * log(p) + (1 - resp) * log(1 - p))
    sum(weights * dev) / sum(weights)
  } else {
    sum(weights * (resp - pred)^2) / sum(weights)
  }
}

# ----------------------------------------------------------------------------
#  predict.meep() -- predictions on NEW data
# ----------------------------------------------------------------------------

#' Predict from a fitted `meep` object
#'
#' Out-of-fold predictions exist only for the training rows (and are stored
#' on the fitted object). For genuinely new data, `predict.meep()` refits
#' each base learner on the **full** training sample -- using the frozen
#' hyperparameters when `tune = "once"` -- then combines the refits with the
#' stored ensemble weights. Full-data refits are cached lazily in an
#' environment field on the object and are not serialized.
#'
#' Predictions on new data can extrapolate; the cross-fitting honesty
#' guarantee applies only to the training-row OOF predictions.
#'
#' @param object A fitted `"meep"` object.
#' @param newdata A matrix or data.frame of new covariates.
#' @param nuisance Which nuisance to predict: `"outcome"`, `"treatment"`,
#'   `"mu0"`, or `"mu1"`.
#' @param ... Unused.
#' @return A numeric vector of length `nrow(newdata)`.
#' @export
predict.meep <- function(object, newdata,
                         nuisance = c("outcome", "treatment",
                                      "mu0", "mu1"),
                         ...) {
  nuisance <- match.arg(nuisance)
  if (missing(newdata) || is.null(newdata))
    stop("meep: predict() requires `newdata`.", call. = FALSE)
  if (!nuisance %in% names(object$oof_matrix))
    stop("meep: nuisance '", nuisance,
         "' was not fitted by this meep() call.", call. = FALSE)

  wts <- object$ensemble_weights[[nuisance]]
  specs <- object$.train$learner_specs
  fam <- if (nuisance == "treatment") object$treatment_family else
    object$family

  # response + training rows for this nuisance
  if (nuisance == "outcome") {
    resp <- object$.train$y; rows <- NULL
  } else if (nuisance == "treatment") {
    resp <- object$.train$d; rows <- NULL
  } else if (nuisance == "mu0") {
    resp <- object$.train$y; rows <- which(object$.train$d == 0)
  } else {
    resp <- object$.train$y; rows <- which(object$.train$d == 1)
  }
  Xtr <- object$.train$X
  if (!is.null(rows)) {
    Xtr <- if (is.data.frame(Xtr)) Xtr[rows, , drop = FALSE] else
      Xtr[rows, , drop = FALSE]
    resp <- resp[rows]
  }

  L <- length(object$learners)
  P <- matrix(NA_real_, nrow(as.data.frame(newdata)), L)
  for (li in seq_len(L)) {
    lname <- object$learners[li]
    if (wts[li] <= 0) next                # learner contributes nothing
    key <- paste0(nuisance, "::", lname)
    fitted_model <- object$.full_fits[[key]]
    if (is.null(fitted_model)) {
      hp <- if (identical(object$tune, "once"))
        (object$frozen_hyperparams[[lname]] %||% list()) else list()
      fitted_model <- tryCatch(
        specs[[li]]$fit(Xtr, resp, fam, object$weights, hp),
        error = function(e) e)
      if (inherits(fitted_model, "error")) next
      assign(key, fitted_model, envir = object$.full_fits)
    }
    pr <- tryCatch(specs[[li]]$predict(fitted_model, newdata, fam),
                   error = function(e) e)
    if (inherits(pr, "error")) next
    P[, li] <- as.numeric(pr)
  }
  .meep_apply_weights(P, wts)
}

# ----------------------------------------------------------------------------
#  print.meep()
# ----------------------------------------------------------------------------

#' @export
print.meep <- function(x, ...) {
  cat("meep(): cross-fitted causal ensemble\n")
  cat("Call: "); print(x$call)
  cat("\n  n = ", x$n, ",  p = ", x$p, "\n", sep = "")
  K <- length(unique(x$folds))
  cat("  folds: ", K, "-fold cross-fitting\n", sep = "")
  cat("  learners: ", paste(x$learners, collapse = ", "), "\n", sep = "")
  cat("  outcome family: ", x$family, sep = "")
  if (!is.null(x$treatment_family))
    cat("  |  treatment family: ", x$treatment_family, sep = "")
  cat("\n")
  cat("  ensemble: ", x$ensemble, "  |  tune: ", x$tune, "\n", sep = "")

  cat("\nEnsemble weights:\n")
  for (nm in names(x$ensemble_weights)) {
    w <- x$ensemble_weights[[nm]]
    cat("  ", format(nm, width = 9), ": ",
        paste(sprintf("%s=%.3f", names(w), w), collapse = "  "),
        "\n", sep = "")
  }

  cat("\nOOF loss (per nuisance, ensemble):\n")
  y_resp <- x$y_resid + x$y_hat_oof          # reconstruct training response
  d_resp <- if (!is.null(x$d_resid)) x$d_resid + x$d_hat_oof else NULL
  for (nm in names(x$oof_matrix)) {
    fam <- if (nm == "treatment") x$treatment_family else x$family
    pred <- switch(nm,
                   outcome   = x$y_hat_oof,
                   treatment = x$d_hat_oof,
                   mu0       = x$mu0_hat_oof,
                   mu1       = x$mu1_hat_oof)
    resp <- if (nm == "treatment") d_resp else y_resp
    ok <- is.finite(pred) & is.finite(resp)
    lv <- if (any(ok)) .meep_loss(resp[ok], pred[ok], fam) else NA_real_
    metric <- if (identical(fam, "binomial")) "deviance" else "mse"
    cat("  ", format(nm, width = 9), ": ", metric, " = ",
        format(lv, digits = 4), "\n", sep = "")
  }

  nfail <- length(x$fold_failures)
  if (nfail > 0L)
    cat("\n  fold failures logged: ", nfail,
        " (see summary() / $fold_failures)\n", sep = "")
  invisible(x)
}

# ----------------------------------------------------------------------------
#  summary.meep()
# ----------------------------------------------------------------------------

#' @export
summary.meep <- function(object, ...) {
  qprobs <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
  yq <- stats::quantile(object$y_hat_oof, qprobs, na.rm = TRUE)
  dq <- if (!is.null(object$d_hat_oof))
    stats::quantile(object$d_hat_oof, qprobs, na.rm = TRUE) else NULL

  overlap <- NULL
  if (!is.null(object$d_hat_oof) &&
      identical(object$treatment_family, "binomial")) {
    ps <- object$d_hat_oof
    overlap <- list(
      min = min(ps, na.rm = TRUE),
      max = max(ps, na.rm = TRUE),
      n_extreme = sum(ps < 0.01 | ps > 0.99, na.rm = TRUE))
  }

  # fold-failure summary table
  ff <- object$fold_failures
  ff_tab <- NULL
  if (length(ff) > 0L) {
    ff_df <- do.call(rbind, lapply(ff, function(z) data.frame(
      nuisance = z$nuisance, fold = z$fold, learner = z$learner,
      reason = z$reason, stringsAsFactors = FALSE)))
    ff_tab <- ff_df
  }

  out <- list(
    call             = object$call,
    n                = object$n,
    p                = object$p,
    K                = length(unique(object$folds)),
    learners         = object$learners,
    family           = object$family,
    treatment_family = object$treatment_family,
    ensemble         = object$ensemble,
    tune             = object$tune,
    ensemble_weights = object$ensemble_weights,
    learner_cv_perf  = object$learner_cv_perf,
    y_quantiles      = yq,
    d_quantiles      = dq,
    overlap          = overlap,
    fold_failures    = ff_tab,
    n_fold_failures  = length(ff)
  )
  class(out) <- "summary.meep"
  out
}

#' @export
print.summary.meep <- function(x, ...) {
  cat("meep() summary\n")
  cat("Call: "); print(x$call)
  cat("\n  n = ", x$n, ",  p = ", x$p, ",  ", x$K, "-fold\n", sep = "")
  cat("  learners: ", paste(x$learners, collapse = ", "), "\n", sep = "")
  cat("  outcome family: ", x$family, sep = "")
  if (!is.null(x$treatment_family))
    cat("  |  treatment family: ", x$treatment_family, sep = "")
  cat("\n  ensemble: ", x$ensemble, "  |  tune: ", x$tune, "\n", sep = "")

  cat("\nPer learner x nuisance OOF performance:\n")
  print(x$learner_cv_perf, row.names = FALSE)

  cat("\nEnsemble weights:\n")
  for (nm in names(x$ensemble_weights)) {
    w <- x$ensemble_weights[[nm]]
    cat("  ", format(nm, width = 9), ": ",
        paste(sprintf("%s=%.3f", names(w), w), collapse = "  "),
        "\n", sep = "")
  }

  cat("\nCross-fitted E[Y|X] quantiles:\n")
  print(round(x$y_quantiles, 4))
  if (!is.null(x$d_quantiles)) {
    cat("\nCross-fitted E[D|X] quantiles:\n")
    print(round(x$d_quantiles, 4))
  }

  if (!is.null(x$overlap)) {
    cat("\nOverlap diagnostic (propensity score):\n")
    cat("  range = [", format(x$overlap$min, digits = 4), ", ",
        format(x$overlap$max, digits = 4), "]\n", sep = "")
    cat("  n outside [0.01, 0.99] = ", x$overlap$n_extreme, sep = "")
    if (x$overlap$n_extreme > 0L)
      cat("  (overlap warning: extreme propensities)")
    cat("\n")
  }

  if (x$n_fold_failures > 0L) {
    cat("\nFold failures (", x$n_fold_failures, " total):\n", sep = "")
    print(x$fold_failures, row.names = FALSE)
  } else {
    cat("\nNo fold failures.\n")
  }
  invisible(x)
}
