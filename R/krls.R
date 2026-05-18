# krls -- Kernel Regularized Least Squares.
#
# File layout:
#   1. krls()              -- S3 generic dispatch.
#   2. krls.formula()      -- formula -> matrix dispatch (mirrors ares.formula).
#   3. krls.default()      -- main fitter, mirrors KRLS::krls (Hainmueller and
#                             Hazlett 2014) numerically.  Returns an S3
#                             object of class c("krls_rr", "krls").
#   4. predict.krls_rr()   -- predictions on new data, optional fit SEs.
#                             Handles formula-side terms / factor expansion.
#   5. print.krls_rr()     -- compact one-screen summary.
#   6. summary.krls_rr()   -- richer summary with quartiles of marginal
#                             effects (mirrors KRLS::summary.krls layout).
#   7. .krls_lambdasearch  -- golden-section LOO lambda selector.
#   8. .krls_fd_binary     -- replace binary-column derivatives with
#                             finite-difference contrasts (mirrors
#                             KRLS::fdskrls).
#
# Conventions:
#   - Functions prefixed with `.` are internal helpers (not exported).
#   - The C++ engine is reached via Rcpp exports: krls_kernel_cpp,
#     krls_kernel_pred_cpp, krls_eig_cpp, krls_solve_cpp,
#     krls_loo_loss_cpp, krls_deriv_cpp, krls_avg_deriv_var_cpp,
#     krls_vsq_cpp (see src/krls.cpp).
#   - Output field names follow KRLS::krls (`coeffs`, `Looe`, `fitted`,
#     `sigma`, `lambda`, `R2`, `derivatives`, `avgderivatives`,
#     `var.avgderivatives`, `vcov.c`, `vcov.fitted`, `binaryindicator`)
#     so downstream code that consumed KRLS output works unchanged.
#
# Naming note: KRLS exports its own krls() with class "krls".  Loading
# both packages masks one or the other depending on attach order.  Use
# `roadrunner::krls()` or `KRLS::krls()` to disambiguate.

#' Kernel Regularized Least Squares
#'
#' Fits a Kernel Regularized Least Squares model with Gaussian kernel,
#' selecting the ridge penalty by leave-one-out cross-validation via the
#' closed-form identity of Hainmueller and Hazlett (2014).  Marginal
#' effects and their variances are computed by default and returned on
#' the original `(X, y)` scale.
#'
#' Three call styles are supported (mirrors `ares()`):
#'
#'   * `krls(X, y, ...)`  -- matrix / numeric interface (back-compatible).
#'   * `krls(y ~ x1 + x2 + ..., data = df, ...)` -- formula interface
#'     with factor expansion and derived terms (`I(x^2)`, `poly(x, 2)`,
#'     etc).
#'   * `krls(df, y, ...)` where `df` is a data frame with numeric +
#'     factor / character columns -- expands categorical columns via
#'     `model.matrix(~ ., df)` (treatment contrasts, intercept dropped).
#'
#' The numerical pipeline mirrors `KRLS::krls()` exactly:
#'
#' 1. `X` and `y` are standardised (column-centred, unit sd).
#' 2. The Gaussian kernel `K_ij = exp(-||x_i - x_j||^2 / sigma)` is built
#'    in parallel C++.
#' 3. `K` is eigendecomposed.  All subsequent solves use the
#'    eigen-basis closed forms, never inverting `K + lambda I` directly.
#' 4. `lambda` is selected by golden-section search on the closed-form
#'    LOO error sum, with the same `(L, U, tol)` bracket as `KRLS::krls`.
#' 5. Marginal effects use the closed-form identity
#'    `dy/dx_k = -(2/sigma) * (X_k * (K c) - K diag(c) X)_k`, computed
#'    without forming the `n x n` distance matrix.
#' 6. Output is unstandardised back to the original `(X, y)` scale.
#'
#' At a fixed `(sigma, lambda)`, fits agree with `KRLS::krls()` to
#' floating-point precision (typically `< 1e-12` on coefficients,
#' fitted values, and marginal effects for `n <= 1000`).  Wall-clock
#' time is roughly `6-10x` faster than `KRLS::krls()` at `n >= 500`
#' when marginal effects and variance estimates are requested.
#'
#' Memory scales as `O(n^2)`: the kernel and its squared eigenvector
#' matrix are both stored.  Expect about `0.4 * n^2 / 1e6` MB of
#' peak working memory (e.g. ~400MB at `n = 1000`, ~10GB at `n = 5000`).
#'
#' @param X A numeric matrix or data frame of predictors (`n x p`).
#'   Constant columns are rejected. Factor / character columns in a
#'   data frame are expanded via `model.matrix(~ ., x)` (treatment
#'   contrasts, intercept dropped). Missing values are handled via
#'   `na.action` (see below).
#' @param y A numeric response vector or single-column matrix.  Constant
#'   `y` is rejected.
#' @param data Used only by the formula method. A data frame containing
#'   the variables referenced by the formula.
#' @param subset Used only by the formula method. An optional integer
#'   or logical vector restricting rows of `data` used for the fit.
#' @param sigma Gaussian-kernel bandwidth.  Default `ncol(X)`.  Must be
#'   a positive scalar.
#' @param lambda Optional ridge penalty.  If `NULL` (default), selected
#'   by golden-section search on the LOO error.
#' @param derivative Logical.  If `TRUE` (default), compute pointwise
#'   marginal effects and their average per variable.  Requires `vcov`.
#' @param binary Logical.  If `TRUE` (default), columns of `X` with
#'   exactly two unique values are treated as binary and their marginal
#'   effects are replaced by predicted-Y first differences (matches
#'   `KRLS::fdskrls`).
#' @param vcov Logical.  If `TRUE` (default), compute the coefficient
#'   covariance and the variance of average marginal effects.
#' @param weights Optional vector of observation weights (length `n`,
#'   strictly positive). Internally normalised to mean 1. Implements
#'   weighted KRLS via a `D K D` transform where `D = diag(sqrt(w))`.
#'   `weights = rep(1, n)` is byte-identical to the unweighted path.
#' @param L,U Optional lower / upper bracket for the lambda search.  If
#'   `NULL`, defaults follow `KRLS::krls()`.
#' @param tol Tolerance for the lambda golden section.  Default
#'   `1e-3 * n`.
#' @param eigtrunc Optional eigenvalue truncation cutoff in `(0, 1]`.
#'   When set, eigenvalues below `eigtrunc * max(d)` are dropped from
#'   the solve.  `NULL` (default) keeps all eigenvalues.
#' @param na.action How to handle missing values in `X`. `"impute"`
#'   (default) replaces NAs with the column median (stored on the fit
#'   and reapplied at `predict()` time). `"omit"` drops rows with any
#'   NA. Missing `y` is always an error.
#' @param trace Integer. `0` is silent, `> 0` enables progress
#'   diagnostics (currently: prints chosen lambda when `> 1`,
#'   golden-section progress when `> 2`). Replaces `print.level`.
#' @param nthreads Integer. Number of threads to use for the C++ kernel
#'   build / decomposition. `0` (default) means use
#'   `RcppParallel::defaultNumThreads()`.
#' @param print.level Deprecated alias for `trace`. Prefer `trace`.
#' @param ... Currently unused (caught for forward compatibility).
#'
#' @return An object of S3 class `c("krls_rr", "krls")` with components
#'   mirroring `KRLS::krls()`: `K`, `coeffs`, `Looe`, `fitted`, `X`,
#'   `y`, `sigma`, `lambda`, `R2`, `derivatives`, `avgderivatives`,
#'   `var.avgderivatives`, `vcov.c`, `vcov.fitted`, `binaryindicator`.
#'   Formula-method fits additionally carry `call`, `terms`, `xlevels`,
#'   `factor_info`, `na.action`, and `na.medians` for downstream
#'   `predict()`, `update()`, and `model.matrix()` support.
#'
#'   `Looe` follows the `KRLS::krls()` scale convention: it is the sum
#'   of squared leave-one-out residuals on the *standardised* `y` scale
#'   multiplied by `sd(y)`, so its units are `[y^2 / sd_y] = [y]`.  This
#'   is preserved for downstream compatibility with code that consumed
#'   `KRLS::krls()` output; it is **not** the LOO MSE in raw-`y`
#'   squared units.
#'
#' @references Hainmueller, J. and C. Hazlett (2014).  "Kernel
#'   Regularized Least Squares: Reducing Misspecification Bias with a
#'   Flexible and Interpretable Machine Learning Approach."  *Political
#'   Analysis* 22(2):143--168.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(rnorm(n * 3), n, 3)
#' colnames(X) <- c("age", "income", "score")
#' y <- sin(X[, 1]) + 0.5 * X[, 2]^2 - 0.3 * X[, 3] +
#'   rnorm(n, sd = 0.2)
#'
#' fit <- krls(X, y)
#' fit
#' fit$avgderivatives           # average marginal effect per variable
#' summary(fit)
#'
#' ## Formula interface with mixed numeric + factor predictors.
#' df <- data.frame(age = X[, 1], income = X[, 2],
#'                  region = factor(sample(c("N", "S"), n, replace = TRUE)),
#'                  y = y)
#' fit_f <- krls(y ~ age + income + region, data = df)
#'
#' ## Predictions on new data with pointwise SEs.
#' Xnew <- matrix(rnorm(20 * 3), 20, 3)
#' colnames(Xnew) <- colnames(X)
#' pr <- predict(fit, Xnew, se.fit = TRUE)
#' head(pr$fit)
#' head(pr$se.fit)
#'
#' @export
krls <- function(X, ...) UseMethod("krls")

#' @rdname krls
#' @export
krls.formula <- function(X, data = NULL, subset = NULL, ..., y = NULL) {
  formula <- X
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()

  # Reject offset() terms: KRLS has no GLM offset concept; absorbing into
  # the design silently would be wrong.
  pre_tt <- stats::terms(formula, data = data)
  if (length(attr(pre_tt, "offset"))) {
    stop("krls: offset() terms are not supported in the formula. ",
         "Subtract the offset from y before calling krls() if needed.",
         call. = FALSE)
  }

  # Build model frame with na.pass so krls.default applies its own
  # na.action. Subset is honoured via the lm()-style match.call() trick
  # so NSE resolves it against the caller's environment.
  mfcall <- match.call(expand.dots = FALSE)
  mfcall$y <- NULL
  mfcall[[1L]] <- quote(stats::model.frame)
  names(mfcall)[names(mfcall) == "X"] <- "formula"
  mfcall$`...` <- NULL
  mfcall$na.action <- quote(stats::na.pass)
  mf <- eval(mfcall, parent.frame())
  yv <- stats::model.response(mf)
  if (is.null(yv))
    stop("krls: response variable is missing from formula/data.")
  mm <- stats::model.matrix(formula, mf)
  has_int <- "(Intercept)" %in% colnames(mm)
  if (has_int) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  out <- krls.default(X = mm, y = as.numeric(yv), ...)
  out$call <- cl
  tt <- stats::terms(formula, data = data)
  out$terms <- tt
  out$xlevels <- stats::.getXlevels(tt, mf)
  out
}

#' @rdname krls
#' @export
krls.default <- function(X, y,
                         sigma = NULL, lambda = NULL,
                         derivative = TRUE, binary = TRUE, vcov = TRUE,
                         weights = NULL,
                         L = NULL, U = NULL, tol = NULL, eigtrunc = NULL,
                         lambda.method = c("loo", "cv"),
                         lambda.grid = NULL,
                         nfold = 0L, ncross = 1L, stratify = TRUE,
                         seed.cv = NULL, cv.1se = FALSE,
                         na.action = c("impute", "omit"),
                         trace = NULL, nthreads = 0L,
                         print.level = NULL, ...) {
  cl <- match.call()
  na.action <- match.arg(na.action)
  lambda.method <- match.arg(lambda.method)

  ## --- print.level / trace harmonisation ---------------------------
  # Phase 8: `trace` is the canonical name (matches ares).  `print.level`
  # is kept as a deprecated alias for back-compat; we warn once if the
  # user supplies it AND it would shadow an explicit `trace`.
  if (!is.null(print.level)) {
    if (is.null(trace)) {
      trace <- print.level
    } else if (!identical(as.integer(trace), as.integer(print.level))) {
      warning("krls: both `trace` and `print.level` supplied; using `trace`.",
              " (`print.level` is deprecated; prefer `trace`.)",
              call. = FALSE)
    } else {
      warning("krls: `print.level` is deprecated; use `trace` instead.",
              call. = FALSE)
    }
  }
  if (is.null(trace)) trace <- 0L
  trace <- as.integer(trace)

  ## --- argument validation -----------------------------------------
  if (is.null(X) || is.null(y)) stop("X and y are required")
  if (is.factor(y) || is.character(y)) {
    stop("y must be numeric (got ",
         if (is.factor(y)) "factor" else "character", ")")
  }
  if (is.factor(X) || is.character(X)) {
    stop("X must be numeric (got ",
         if (is.factor(X)) "factor" else "character", ")")
  }

  ## --- data.frame factor / character expansion ---------------------
  factor_info <- NULL
  pre_na_medians <- NULL
  if (is.data.frame(X)) {
    is_cat <- vapply(X, function(z) is.factor(z) || is.character(z),
                     logical(1L))
    is_num <- vapply(X, is.numeric, logical(1L))

    # Numeric-column NA handling BEFORE model.matrix (which would drop
    # rows under its default na.action). Mirrors ares.default.
    if (any(is_num)) {
      num_names <- names(X)[is_num]
      df_na_mask <- vapply(num_names,
                           function(jn) anyNA(X[[jn]]), logical(1L))
      if (any(df_na_mask)) {
        rows_aff <- sum(rowSums(is.na(
          as.matrix(X[, num_names, drop = FALSE]))) > 0L)
        cols_aff <- sum(df_na_mask)
        if (na.action == "impute") {
          pre_na_medians <- vapply(
            num_names,
            function(jn) stats::median(X[[jn]], na.rm = TRUE),
            numeric(1L))
          names(pre_na_medians) <- num_names
          if (any(!is.finite(pre_na_medians)))
            stop("krls: na.action='impute' but at least one numeric column",
                 " is entirely NA; cannot compute median.",
                 " Drop the column or pass na.action='omit'.")
          for (jn in num_names) {
            na_j <- is.na(X[[jn]])
            if (any(na_j)) X[[jn]][na_j] <- pre_na_medians[[jn]]
          }
          warning("krls: median-imputed ", sum(df_na_mask),
                  " numeric column(s) across ", rows_aff,
                  " row(s). Column medians are stored on the fit.",
                  call. = FALSE)
        } else {
          keep <- stats::complete.cases(X[, num_names, drop = FALSE])
          dropped <- sum(!keep)
          warning("krls: dropped ", dropped,
                  " incomplete row(s) across ", cols_aff,
                  " numeric column(s). Pass na.action='impute' to",
                  " median-impute instead.", call. = FALSE)
          X <- X[keep, , drop = FALSE]
          y <- y[keep]
          if (!is.null(weights)) weights <- weights[keep]
          if (length(y) < 3L)
            stop("krls: na.action='omit' left fewer than 3 rows; aborting.")
        }
      }
    }

    if (any(is_cat)) {
      for (j in which(is_cat)) {
        if (is.character(X[[j]])) X[[j]] <- factor(X[[j]])
      }
      xlevels <- lapply(X[is_cat], levels)
      xmm <- stats::model.matrix(~ ., data = X)
      if ("(Intercept)" %in% colnames(xmm))
        xmm <- xmm[, colnames(xmm) != "(Intercept)", drop = FALSE]
      factor_info <- list(
        xlevels       = xlevels,
        orig_names    = names(X),
        is_cat        = is_cat,
        is_num        = is_num,
        expanded_names = colnames(xmm),
        num_medians   = pre_na_medians
      )
      X <- xmm
    } else {
      nm <- names(X)
      X <- as.matrix(X)
      storage.mode(X) <- "double"
      colnames(X) <- nm
    }
  }

  X <- as.matrix(X)
  y <- as.matrix(y)
  if (ncol(y) != 1L) {
    stop("y must be a vector or single-column matrix (got ", ncol(y),
         " columns)")
  }
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  if (!is.numeric(X)) stop("X must be numeric")
  if (!is.numeric(y)) stop("y must be numeric")
  if (anyNA(y))      stop("krls: y contains missing data")

  ## --- na.action on numeric matrix ---------------------------------
  na_medians <- NULL
  if (any(is.nan(X)) || any(is.infinite(X)))
    stop("krls: X contains NaN or +/-Inf values; aborting fit.",
         " (NA is handled via `na.action`; NaN / Inf are not.)")
  na_x <- is.na(X)
  if (any(na_x)) {
    rows_aff <- sum(rowSums(na_x) > 0L)
    cols_aff <- sum(colSums(na_x) > 0L)
    if (na.action == "impute") {
      na_medians <- vapply(seq_len(ncol(X)),
                           function(j) stats::median(X[, j], na.rm = TRUE),
                           numeric(1L))
      names(na_medians) <- colnames(X)
      if (any(!is.finite(na_medians)))
        stop("krls: na.action='impute' but at least one X column is",
             " entirely NA; cannot compute median.")
      for (j in seq_len(ncol(X))) {
        na_j <- na_x[, j]
        if (any(na_j)) X[na_j, j] <- na_medians[j]
      }
      warning("krls: median-imputed ", sum(na_x), " NA value(s) across ",
              rows_aff, " row(s) and ", cols_aff, " column(s).",
              " Column medians stored on fit.", call. = FALSE)
    } else {
      keep <- rowSums(na_x) == 0L
      dropped <- sum(!keep)
      warning("krls: dropped ", dropped, " incomplete row(s) across ",
              cols_aff, " column(s).", call. = FALSE)
      X <- X[keep, , drop = FALSE]
      y <- y[keep, , drop = FALSE]
      if (!is.null(weights)) weights <- weights[keep]
      if (length(y) < 3L)
        stop("krls: na.action='omit' left fewer than 3 rows; aborting.")
    }
  }

  if (var(as.vector(y)) == 0) stop("y is a constant (does not vary)")
  stopifnot(is.logical(derivative), is.logical(vcov), is.logical(binary))
  if (derivative && !vcov) {
    stop("derivative = TRUE requires vcov = TRUE")
  }
  n <- nrow(X); d <- ncol(X)
  if (n != nrow(y)) stop("nrow(X) not equal to length of y")

  ## --- weights validation (Phase 2) -------------------------------
  w_norm <- NULL
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("krls: weights must be numeric.")
    if (length(weights) != n)
      stop("krls: length(weights) (", length(weights),
           ") must equal nrow(X) (", n, ").")
    if (any(!is.finite(weights)))
      stop("krls: weights contain NA / NaN / Inf; all weights must be finite.")
    if (any(weights <= 0))
      stop("krls: weights must be strictly positive.")
    wm <- mean(weights)
    if (wm <= 0)
      stop("krls: weights mean is zero; cannot normalise.")
    w_norm <- as.numeric(weights) / wm
  }

  if (!is.null(eigtrunc)) {
    stopifnot(is.numeric(eigtrunc), length(eigtrunc) == 1)
    if (eigtrunc < 0 || eigtrunc > 1) {
      stop("eigtrunc must be in [0, 1]")
    }
    if (eigtrunc == 0) {
      eigtrunc <- NULL
      warning("eigtrunc = 0 is equivalent to NULL; ignoring")
    }
  }
  if (is.null(sigma)) {
    sigma <- d
  } else {
    stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(d))

  ## --- nthreads (Phase 8) -----------------------------------------
  nthreads <- as.integer(nthreads)
  if (is.na(nthreads) || nthreads < 0L) nthreads <- 0L
  nthreads_eff <- if (nthreads <= 0L)
    RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  ## --- standardise -------------------------------------------------
  X.init     <- X
  X.init.sd  <- apply(X.init, 2L, sd)
  if (any(X.init.sd == 0)) {
    stop("at least one column in X is a constant, please remove the constant(s)")
  }
  y.init      <- y
  y.init.sd   <- apply(y.init, 2L, sd)
  y.init.mean <- mean(y.init)
  Xs <- scale(X.init, center = TRUE, scale = X.init.sd)
  ys <- scale(y.init, center = y.init.mean, scale = y.init.sd)
  ## drop attributes for clean C++ handoff
  Xs <- matrix(Xs, n, d, dimnames = list(NULL, colnames(X)))
  ys <- as.numeric(ys)

  ## --- weighted-LS transform (Phase 2) ------------------------------
  # Solve a *weighted* KRLS by applying D = diag(sqrt(w)) on both sides
  # of K. The closed-form fit on D K D with target D ys recovers the
  # weighted normal equations
  #   (K W K + lambda K) c = K W y   <=>   (K + lambda K^-1) c = W y
  # in the standardised space. Coefficients are returned in the
  # original (un-D-twisted) basis: c_final = D c_solve so that
  # K %*% c_final == K_w %*% c_solve and yhat is on the un-weighted
  # scale (n x 1 prediction on the training X). When weights = rep(1, n)
  # this collapses byte-identically to the unweighted path (sqw is the
  # vector of ones; D is the identity).
  if (!is.null(w_norm)) {
    sqw <- sqrt(w_norm)
  } else {
    sqw <- NULL
  }

  ## --- kernel + eigendecomposition ---------------------------------
  K <- krls_kernel_cpp(Xs, sigma)
  if (!is.null(sqw)) {
    # Twist the kernel: K_w = D K D
    K_w <- K * tcrossprod(sqw)
    eo <- krls_eig_cpp(K_w)
    dvals <- as.numeric(eo$values)
    V     <- eo$vectors
    Vsq   <- krls_vsq_cpp(V)
    ys_w  <- sqw * ys
    Vty   <- as.numeric(crossprod(V, ys_w))
  } else {
    eo <- krls_eig_cpp(K)
    dvals <- as.numeric(eo$values)
    V     <- eo$vectors
    Vsq   <- krls_vsq_cpp(V)
    Vty   <- as.numeric(crossprod(V, ys))
  }

  ## --- eigentruncation handling (optional) -------------------------
  if (!is.null(eigtrunc)) {
    lastkeeper <- max(which(dvals >= eigtrunc * dvals[1L]))
    lastkeeper <- max(1L, lastkeeper)
    if (lastkeeper < length(dvals)) {
      dvals <- dvals[seq_len(lastkeeper)]
      V     <- V[, seq_len(lastkeeper), drop = FALSE]
      Vsq   <- Vsq[, seq_len(lastkeeper), drop = FALSE]
      Vty   <- Vty[seq_len(lastkeeper)]
    }
  }

  ## --- lambda selection --------------------------------------------
  # Phase 3: lambda.method = c('loo', 'cv'). LOO is the back-compat
  # default and uses the closed-form Hainmueller-Hazlett identity. CV
  # does K-fold (or repeated K-fold) over an explicit lambda.grid and
  # picks argmin held-out MSE (or 1-SE), optionally seeded for
  # reproducibility.
  cv_diag <- NULL
  if (is.null(lambda)) {
    if (lambda.method == "loo") {
      lambda <- .krls_lambdasearch(dvals = dvals, V = V, Vsq = Vsq,
                                   Vty = Vty, n_y = nrow(y),
                                   L = L, U = U, tol = tol,
                                   noisy = isTRUE(trace > 2))
      if (trace > 1) {
        cat("Lambda that minimizes Loo-Loss is:", round(lambda, 5), "\n")
      }
    } else {
      # CV path. Default nfold = 5 if user asked lambda.method='cv' but
      # didn't specify nfold. nfold = nrow(X) reproduces LOO via CV.
      nfold_int <- as.integer(nfold)
      if (is.na(nfold_int) || nfold_int <= 0L) nfold_int <- 5L
      if (nfold_int > nrow(X))
        stop("krls: nfold (", nfold_int, ") exceeds nrow(X) (",
             nrow(X), ").")
      ncross_int <- as.integer(ncross)
      if (is.na(ncross_int) || ncross_int < 1L) ncross_int <- 1L

      # Default lambda.grid: 30 log-spaced points within the standard
      # (L, U) bracket from the closed-form search (which is itself
      # derived from the eigenvalues so the grid is data-adaptive).
      if (is.null(lambda.grid)) {
        # Reuse the same L/U bracket the golden search would use.
        U_grid <- if (is.null(U)) {
          Uloc <- nrow(y)
          while (sum(dvals / (dvals + Uloc)) < 1) Uloc <- Uloc - 1
          Uloc
        } else U
        L_grid <- if (is.null(L)) {
          q  <- which.min(abs(dvals - (max(dvals) / 1000)))
          Lloc  <- .Machine$double.eps
          while (sum(dvals / (dvals + Lloc)) > q) Lloc <- Lloc + 0.05
          Lloc
        } else max(L, .Machine$double.eps)
        if (U_grid <= L_grid) {
          U_grid <- max(L_grid * 10, L_grid + 1)
        }
        lambda.grid <- exp(seq(log(max(L_grid, .Machine$double.eps)),
                               log(U_grid), length.out = 30L))
      } else {
        stopifnot(is.numeric(lambda.grid),
                  length(lambda.grid) >= 1L,
                  all(lambda.grid > 0))
        lambda.grid <- sort(unique(as.numeric(lambda.grid)))
      }

      # RNG protection so CV only consumes the user's RNG when
      # seed.cv = NULL.
      if (!is.null(seed.cv)) {
        if (exists(".Random.seed", envir = globalenv(),
                   inherits = FALSE)) {
          old_seed <- get(".Random.seed", envir = globalenv(),
                          inherits = FALSE)
          on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                  add = TRUE)
        } else {
          on.exit(rm(list = ".Random.seed", envir = globalenv()),
                  add = TRUE)
        }
        set.seed(as.integer(seed.cv))
      }

      n_cv <- nrow(X)
      # Stratification for continuous y: use quantile bins (10 bins or
      # fewer if n < 50).  For a small-y, fall back to plain shuffle.
      build_folds <- function(n_pts, k, strat) {
        if (isTRUE(strat) && n_pts >= 20L) {
          nb <- min(10L, max(2L, floor(n_pts / 5)))
          bins <- cut(as.numeric(y), breaks = nb, include.lowest = TRUE,
                      labels = FALSE)
          # Distribute each bin's indices round-robin across folds.
          fid <- integer(n_pts)
          for (b in seq_len(nb)) {
            idx_b <- which(bins == b)
            if (length(idx_b) == 0L) next
            idx_b <- sample(idx_b)  # shuffle within bin
            fid[idx_b] <- rep_len(seq_len(k), length(idx_b))
          }
          # Any unassigned (NA bin) get round-robin too.
          if (any(fid == 0L)) {
            miss <- which(fid == 0L)
            fid[miss] <- rep_len(seq_len(k), length(miss))
          }
          fid
        } else {
          rep_len(seq_len(k), n_pts)[sample.int(n_pts)]
        }
      }

      mse_mat <- matrix(NA_real_, nrow = length(lambda.grid),
                        ncol = nfold_int * ncross_int)
      col_idx <- 1L
      for (cc in seq_len(ncross_int)) {
        fid <- build_folds(n_cv, nfold_int, stratify)
        for (k in seq_len(nfold_int)) {
          test_i  <- which(fid == k)
          train_i <- which(fid != k)
          if (length(test_i) == 0L || length(train_i) < 3L) next
          Xs_tr <- Xs[train_i, , drop = FALSE]
          ys_tr <- ys[train_i]
          Xs_te <- Xs[test_i,  , drop = FALSE]
          ys_te <- ys[test_i]
          # Fold-local kernel + eigendecomp (the shared piece across
          # all lambda candidates per fold).
          K_tr <- krls_kernel_cpp(Xs_tr, sigma)
          K_te <- krls_kernel_pred_cpp(Xs_te, Xs_tr, sigma)
          if (!is.null(w_norm)) {
            sqw_tr <- sqrt(w_norm[train_i])
            K_tr_use <- K_tr * tcrossprod(sqw_tr)
            ys_tr_use <- sqw_tr * ys_tr
          } else {
            sqw_tr <- NULL
            K_tr_use <- K_tr
            ys_tr_use <- ys_tr
          }
          eo_k <- krls_eig_cpp(K_tr_use)
          dvals_k <- as.numeric(eo_k$values)
          V_k     <- eo_k$vectors
          Vsq_k   <- krls_vsq_cpp(V_k)
          Vty_k   <- as.numeric(crossprod(V_k, ys_tr_use))
          for (li in seq_along(lambda.grid)) {
            lam_l <- lambda.grid[li]
            sol_k <- krls_solve_cpp(dvals_k, V_k, Vsq_k, Vty_k, lam_l)
            c_solve <- as.numeric(sol_k$coeffs)
            c_fold  <- if (!is.null(sqw_tr)) sqw_tr * c_solve else c_solve
            yhat_te <- as.numeric(K_te %*% c_fold)
            err_v   <- (ys_te - yhat_te)
            if (!is.null(w_norm)) {
              w_te <- w_norm[test_i]
              mse_mat[li, col_idx] <-
                sum(w_te * err_v^2) / max(sum(w_te), .Machine$double.eps)
            } else {
              mse_mat[li, col_idx] <- mean(err_v^2)
            }
          }
          col_idx <- col_idx + 1L
        }
      }
      # Average across folds / repeats. Use min().
      mean_mse <- rowMeans(mse_mat, na.rm = TRUE)
      if (all(is.na(mean_mse)))
        stop("krls: lambda.method='cv' produced no usable folds; ",
             "check nfold / stratify.")
      best_i <- which.min(mean_mse)
      lambda <- lambda.grid[best_i]
      if (isTRUE(cv.1se)) {
        # 1-SE rule: pick the largest lambda whose mean is within 1 SE
        # of the best. SE = sd across folds / sqrt(n_folds).
        sd_mse <- apply(mse_mat, 1L, function(z) sd(z, na.rm = TRUE))
        nfok <- rowSums(!is.na(mse_mat))
        se_mse <- sd_mse / sqrt(pmax(nfok, 1L))
        thresh <- mean_mse[best_i] + se_mse[best_i]
        cand <- which(mean_mse <= thresh)
        if (length(cand)) lambda <- max(lambda.grid[cand])
      }
      cv_diag <- list(
        method      = "cv",
        nfold       = nfold_int,
        ncross      = ncross_int,
        stratify    = isTRUE(stratify),
        seed.cv     = seed.cv,
        lambda.grid = lambda.grid,
        mean_mse    = mean_mse,
        mse_per_fold = mse_mat,
        best_idx    = best_i,
        cv.1se      = isTRUE(cv.1se)
      )
      if (trace > 1) {
        cat("CV lambda (", nfold_int, "-fold x ", ncross_int, " cross):",
            round(lambda, 6), "\n", sep = "")
      }
    }
  } else {
    stopifnot(is.numeric(lambda), length(lambda) == 1, lambda > 0)
  }

  ## --- solve for c, LOOe, diag(G^{-1}) -----------------------------
  sol     <- krls_solve_cpp(dvals, V, Vsq, Vty, lambda)
  coeffs_solve <- as.numeric(sol$coeffs)
  Le      <- sol$Le
  if (!is.null(sqw)) {
    # Un-twist coefficients: c_final = D c_solve so K c_final has the
    # right meaning on the original (un-D-scaled) input.
    coeffs <- sqw * coeffs_solve
    yfit_s <- as.numeric(K %*% coeffs)
  } else {
    coeffs <- coeffs_solve
    yfit_s  <- as.numeric(K %*% coeffs)
  }

  ## --- coefficient covariance --------------------------------------
  vcov.c <- NULL
  vcov.fitted <- NULL
  if (!is.null(sqw)) {
    # Weighted residual variance: sum w * (ys - yfit_s)^2 / sum(w).
    sigma2 <- as.numeric(sum(w_norm * (ys - yfit_s)^2) / sum(w_norm))
  } else {
    sigma2 <- as.numeric((1 / n) * crossprod(ys - yfit_s))
  }
  if (isTRUE(vcov)) {
    w <- sigma2 / (dvals + lambda)^2
    vcovmatc_solve <- V %*% (t(V) * w)
    if (!is.null(sqw)) {
      # Un-twist: vcov(c_final) = D %*% vcov(c_solve) %*% D
      vcovmatc <- vcovmatc_solve * tcrossprod(sqw)
    } else {
      vcovmatc <- vcovmatc_solve
    }
    vcovmatyhat <- crossprod(K, vcovmatc %*% K)
  }

  ## --- derivatives --------------------------------------------------
  derivmat <- avgderiv <- varavgderivmat <- NULL
  if (isTRUE(derivative)) {
    derivmat_s <- krls_deriv_cpp(Xs, K, coeffs, sigma)
    if (!is.null(sqw)) {
      avgderiv_s <- matrix(
        colSums(derivmat_s * w_norm) / sum(w_norm),
        nrow = 1L)
    } else {
      avgderiv_s <- matrix(colMeans(derivmat_s), nrow = 1L)
    }
    varavg_s   <- krls_avg_deriv_var_cpp(Xs, K, V, dvals, sigma, lambda, sigma2)
    colnames(derivmat_s)  <- colnames(X)
    colnames(avgderiv_s)  <- colnames(X)
    scale_vec <- as.numeric(y.init.sd) / X.init.sd
    derivmat  <- sweep(derivmat_s, 2L, scale_vec, "*")
    if (!is.null(sqw)) {
      avgderiv  <- matrix(
        colSums(derivmat * w_norm) / sum(w_norm),
        nrow = 1L,
        dimnames = list(NULL, colnames(X)))
    } else {
      avgderiv  <- matrix(colMeans(derivmat), nrow = 1L,
                          dimnames = list(NULL, colnames(X)))
    }
    varavgderivmat <- matrix(scale_vec^2 * as.numeric(varavg_s),
                             nrow = 1L,
                             dimnames = list(NULL, colnames(X)))
  }

  ## --- unstandardise fitted, vcov ----------------------------------
  yfitted <- yfit_s * as.numeric(y.init.sd) + y.init.mean
  if (isTRUE(vcov)) {
    vcov.c      <- (as.numeric(y.init.sd)^2) * vcovmatc
    vcov.fitted <- (as.numeric(y.init.sd)^2) * vcovmatyhat
  }

  Looe <- Le * as.numeric(y.init.sd)
  R2   <- 1 - (var(as.numeric(y.init - yfitted)) / (as.numeric(y.init.sd)^2))

  binaryindicator <- matrix(FALSE, 1L, d,
                            dimnames = list(NULL, colnames(X)))

  out <- list(K = K, coeffs = coeffs, Looe = Looe, fitted = yfitted,
              X = X.init, y = y.init,
              sigma = sigma, lambda = lambda, R2 = R2,
              derivatives = derivmat,
              avgderivatives = avgderiv,
              var.avgderivatives = varavgderivmat,
              vcov.c = vcov.c, vcov.fitted = vcov.fitted,
              binaryindicator = binaryindicator)

  # Bookkeeping for formula method / predict() / update()
  out$call <- cl
  out$na.action <- na.action
  out$na.medians <- na_medians
  out$factor_info <- factor_info
  out$weights <- if (!is.null(w_norm)) w_norm else NULL
  if (!is.null(cv_diag)) out$cv <- cv_diag

  class(out) <- c("krls_rr", "krls")

  if (isTRUE(derivative) && isTRUE(binary)) {
    out <- .krls_fd_binary(out)
  }

  if (trace > 0 && isTRUE(derivative)) {
    av <- setNames(as.vector(out$avgderivatives),
                   colnames(out$avgderivatives))
    cat("\n Average Marginal Effects:\n \n"); print(av)
    cat("\n Quartiles of Marginal Effects:\n \n")
    print(apply(out$derivatives, 2L, quantile,
                probs = c(0.25, 0.5, 0.75)))
  }
  out
}

#' @rdname krls
#' @param object A fitted `"krls"` object.
#' @param newdata A numeric matrix or data frame with the same columns
#'   as the training `X` (or matching the training formula).  `NULL`
#'   returns `object$fitted`.
#' @param se.fit Logical.  If `TRUE`, return pointwise standard errors
#'   of the predictions.  Requires the fit was created with
#'   `vcov = TRUE`.
#' @export
predict.krls_rr <- function(object, newdata = NULL, se.fit = FALSE, ...) {
  if (!inherits(object, "krls")) {
    stop("object is not of class 'krls'")
  }
  if (isTRUE(se.fit) && is.null(object$vcov.c)) {
    stop("refit with krls(..., vcov = TRUE) to compute standard errors")
  }

  if (is.null(newdata)) {
    return(list(fit = matrix(object$fitted, ncol = 1L),
                se.fit = NULL, vcov.fit = NULL,
                newdata = NULL, newdataK = NULL))
  }

  # Build the design matrix to score against.
  xnew <- .krls_build_design(object, newdata)

  storage.mode(xnew) <- "double"
  if (ncol(object$X) != ncol(xnew)) {
    stop("ncol(newdata) (", ncol(xnew),
         ") differs from ncol(X) from fitted krls object (",
         ncol(object$X), ")")
  }

  # Re-impute NAs in newdata using stored medians (matrix path only;
  # the data-frame / formula path imputed inside .krls_build_design).
  if (any(is.na(xnew))) {
    if (!is.null(object$na.medians)) {
      for (j in seq_len(ncol(xnew))) {
        na_j <- is.na(xnew[, j])
        if (any(na_j)) xnew[na_j, j] <- object$na.medians[j]
      }
      warning("krls: median-imputed ", sum(is.na(xnew)),
              " missing value(s) in newdata using stored training medians.",
              call. = FALSE)
    } else {
      stop("krls: newdata has NA values but the fit has no stored medians",
           " (training used na.action='omit'). Drop or impute the rows",
           " yourself before calling predict().")
    }
  }

  # Column-name handling for the bare-matrix path (only applies when
  # no factor_info / terms drove the rebuild above).
  if (is.null(object$factor_info) && is.null(object$terms)) {
    cn_train <- colnames(object$X)
    cn_new   <- colnames(xnew)
    if (!is.null(cn_train) && !is.null(cn_new)) {
      if (!setequal(cn_train, cn_new)) {
        stop("colnames(newdata) do not match colnames(X) from fit: ",
             "training = c(", paste(shQuote(cn_train), collapse = ", "),
             "), newdata = c(", paste(shQuote(cn_new), collapse = ", "),
             ")")
      }
      if (!identical(cn_train, cn_new)) {
        xnew <- xnew[, cn_train, drop = FALSE]
      }
    }
  }

  if (nrow(object$X) < 2L) {
    stop("object$X has fewer than 2 rows; sd cannot be recomputed")
  }
  Xmeans <- colMeans(object$X)
  Xsd    <- apply(object$X, 2L, sd)
  if (any(!is.finite(Xsd)) || any(Xsd == 0)) {
    stop("at least one stored X column has zero or non-finite sd; ",
         "object$X may have been mutated after fit")
  }
  Xs     <- scale(object$X, center = Xmeans, scale = Xsd)
  Xs     <- matrix(Xs, nrow(object$X), ncol(object$X))
  Xn     <- scale(xnew, center = Xmeans, scale = Xsd)
  Xn     <- matrix(Xn, nrow(xnew), ncol(xnew))

  Knew   <- krls_kernel_pred_cpp(Xn, Xs, object$sigma)
  yhat   <- Knew %*% object$coeffs

  vcov.fit <- se.fit.out <- NULL
  if (isTRUE(se.fit)) {
    vcov.c.raw  <- object$vcov.c * as.vector(1 / var(as.vector(object$y)))
    vcov.fitted <- tcrossprod(Knew %*% vcov.c.raw, Knew)
    sd_y2       <- as.numeric(apply(object$y, 2L, sd))^2
    vcov.fit    <- sd_y2 * vcov.fitted
    se.fit.out  <- matrix(sqrt(diag(vcov.fit)), ncol = 1L)
  }
  yhat <- (yhat * as.numeric(apply(object$y, 2L, sd))) + mean(object$y)
  list(fit = yhat, se.fit = se.fit.out, vcov.fit = vcov.fit,
       newdata = Xn, newdataK = Knew)
}

## -----------------------------------------------------------------
## Helpers (not exported).
## -----------------------------------------------------------------

# Build a numeric design matrix from `newdata`, honouring the fit's
# stored factor / terms / column metadata. Ports BUG-004 (OOV factor),
# BUG-010 (factor NA), and BUG-012 (derived terms) from predict.ares.
.krls_build_design <- function(object, newdata) {
  if (is.data.frame(newdata)) {
    if (!is.null(object$factor_info)) {
      fi <- object$factor_info
      if (!all(fi$orig_names %in% colnames(newdata)))
        stop("krls: newdata is missing columns: ",
             paste(setdiff(fi$orig_names, colnames(newdata)),
                   collapse = ", "))
      newdata <- newdata[, fi$orig_names, drop = FALSE]
      # Numeric NA fill with stored training medians.
      if (!is.null(fi$num_medians)) {
        for (jn in names(fi$num_medians)) {
          if (jn %in% names(newdata)) {
            col <- newdata[[jn]]
            if (is.numeric(col) && anyNA(col)) {
              col[is.na(col)] <- fi$num_medians[[jn]]
              newdata[[jn]] <- col
            }
          }
        }
      }
      # Factor NA detection (BUG-010 style).
      fna_report <- character(0)
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
            sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                    jname, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("krls: NA value(s) in factor/character newdata column(s); ",
             "drop the row(s) or impute first.\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      # OOV detection (BUG-004 style).
      oov_report <- list()
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if (is.character(col) || is.factor(col)) {
          col_chr <- as.character(col)
          tr_lev <- fi$xlevels[[jname]]
          bad_rows <- which(!is.na(col_chr) & !(col_chr %in% tr_lev))
          if (length(bad_rows)) {
            bad_lev <- unique(col_chr[bad_rows])
            oov_report[[jname]] <- list(rows = bad_rows, levels = bad_lev)
          }
          newdata[[jname]] <- factor(col, levels = tr_lev)
        }
      }
      if (length(oov_report)) {
        msgs <- vapply(names(oov_report), function(jn) {
          rep <- oov_report[[jn]]
          sprintf("  column %s: %d row(s) (e.g. row %d) with level(s) not seen at fit time: %s",
                  jn, length(rep$rows), rep$rows[1],
                  paste(utils::head(rep$levels, 5), collapse = ", "))
        }, character(1L))
        stop("krls: out-of-vocabulary factor level(s) in newdata.\n",
             paste(msgs, collapse = "\n"), call. = FALSE)
      }
      xnew <- stats::model.matrix(~ ., data = newdata)
      if ("(Intercept)" %in% colnames(xnew))
        xnew <- xnew[, colnames(xnew) != "(Intercept)", drop = FALSE]
      if (!identical(colnames(xnew), fi$expanded_names)) {
        missing_cols <- setdiff(fi$expanded_names, colnames(xnew))
        if (length(missing_cols))
          stop("krls: newdata expansion is missing columns: ",
               paste(missing_cols, collapse = ", "),
               " -- likely an out-of-vocabulary factor level.")
        xnew <- xnew[, fi$expanded_names, drop = FALSE]
      }
      return(as.matrix(xnew))
    }
    if (!is.null(object$terms)) {
      tt <- stats::delete.response(object$terms)
      term_vars <- all.vars(tt)
      fna_report <- character(0)
      for (jn in intersect(term_vars, colnames(newdata))) {
        col <- newdata[[jn]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
            sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                    jn, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("krls: NA value(s) in factor/character newdata column(s).\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      xlev <- if (!is.null(object$xlevels)) object$xlevels else NULL
      mm <- stats::model.matrix(tt, data = newdata, xlev = xlev)
      if ("(Intercept)" %in% colnames(mm))
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      # OOV factor: model.matrix with xlev silently coerces unseen levels
      # to NA in the dummy columns, then na.action=na.omit drops rows.
      # Detect mismatch in row count and stop.
      if (nrow(mm) != nrow(newdata)) {
        stop("krls: newdata expansion dropped rows (likely NA / OOV factor",
             " levels in formula terms). nrow(newdata)=", nrow(newdata),
             ", nrow(mm)=", nrow(mm), ".", call. = FALSE)
      }
      tr_names <- colnames(object$X)
      if (!is.null(tr_names) && !identical(colnames(mm), tr_names)) {
        missing_cols <- setdiff(tr_names, colnames(mm))
        if (length(missing_cols))
          stop("krls: newdata model-matrix expansion is missing columns: ",
               paste(missing_cols, collapse = ", "), call. = FALSE)
        mm <- mm[, tr_names, drop = FALSE]
      }
      return(as.matrix(mm))
    }
    # Plain df, no factor_info, no terms: coerce to matrix.
    return(as.matrix(newdata))
  }
  as.matrix(newdata)
}

.krls_lambdasearch <- function(dvals, V, Vsq, Vty, n_y,
                               L = NULL, U = NULL, tol = NULL,
                               noisy = FALSE) {
  n <- n_y
  if (is.null(tol)) {
    tol <- 1e-3 * n
  } else {
    stopifnot(is.numeric(tol), length(tol) == 1, tol > 0)
  }
  if (is.null(U)) {
    U <- n
    while (sum(dvals / (dvals + U)) < 1) U <- U - 1
  } else {
    stopifnot(is.numeric(U), length(U) == 1, U > 0)
  }
  if (is.null(L)) {
    q  <- which.min(abs(dvals - (max(dvals) / 1000)))
    L  <- .Machine$double.eps
    while (sum(dvals / (dvals + L)) > q) L <- L + 0.05
  } else {
    stopifnot(is.numeric(L), length(L) == 1, L >= 0)
  }
  gr <- 0.381966
  X1 <- L + gr * (U - L)
  X2 <- U - gr * (U - L)
  S1 <- krls_loo_loss_cpp(dvals, V, Vsq, Vty, X1)
  S2 <- krls_loo_loss_cpp(dvals, V, Vsq, Vty, X2)
  if (noisy) {
    cat("L:", L, "X1:", X1, "X2:", X2, "U:", U, "S1:", S1, "S2:", S2, "\n")
  }
  while (abs(S1 - S2) > tol) {
    if (S1 < S2) {
      U  <- X2; X2 <- X1
      X1 <- L + gr * (U - L)
      S2 <- S1
      S1 <- krls_loo_loss_cpp(dvals, V, Vsq, Vty, X1)
    } else {
      L  <- X1; X1 <- X2
      X2 <- U - gr * (U - L)
      S1 <- S2
      S2 <- krls_loo_loss_cpp(dvals, V, Vsq, Vty, X2)
    }
    if (noisy) {
      cat("L:", L, "X1:", X1, "X2:", X2, "U:", U,
          "S1:", S1, "S2:", S2, "\n")
    }
  }
  if (S1 < S2) X1 else X2
}

.krls_fd_binary <- function(object) {
  d <- ncol(object$X); n <- nrow(object$X)
  lu <- function(x) length(unique(x))
  binidx <- which(apply(object$X, 2L, lu) == 2L)
  if (length(binidx) == 0L) return(object)
  est  <- se <- matrix(NA_real_, nrow = 1L, ncol = length(binidx))
  diffs_store <- matrix(NA_real_, nrow = n, ncol = length(binidx))
  for (i in seq_along(binidx)) {
    X1 <- X0 <- object$X
    X1[, binidx[i]] <- max(X1[, binidx[i]])
    X0[, binidx[i]] <- min(X0[, binidx[i]])
    Xall <- rbind(X1, X0)
    h <- matrix(rep(c(1 / n, -(1 / n)), each = n), ncol = 1L)
    pout <- predict.krls_rr(object, newdata = Xall, se.fit = TRUE)
    est[1L, i] <- as.numeric(t(h) %*% pout$fit)
    se[1L, i]  <- as.numeric(sqrt(t(h) %*% pout$vcov.fit %*% h)) * sqrt(2)
    diffs_store[, i] <- pout$fit[1:n] - pout$fit[(n + 1):(2 * n)]
  }
  object$derivatives[, binidx]        <- diffs_store
  object$avgderivatives[, binidx]     <- est
  object$var.avgderivatives[, binidx] <- se^2
  object$binaryindicator[, binidx]    <- TRUE
  object
}

#' @export
print.krls_rr <- function(x, ...) {
  cat("Kernel Regularized Least Squares (KRLS)\n")
  cat("  n =", nrow(x$X), "  p =", ncol(x$X), "\n")
  cat("  sigma =", signif(x$sigma, 4),
      "  lambda =", signif(x$lambda, 4),
      "  R^2 =", signif(x$R2, 4), "\n")
  if (!is.null(x$weights)) {
    cat("  weighted: yes (mean(w)=1 normalisation)\n")
  }
  if (!is.null(x$avgderivatives)) {
    cat("\nAverage Marginal Effects:\n")
    print(setNames(as.vector(x$avgderivatives),
                   colnames(x$avgderivatives)))
  }
  invisible(x)
}

#' @export
summary.krls_rr <- function(object, ...) {
  out <- list()
  out$call_info <- c(n = nrow(object$X), p = ncol(object$X),
                     sigma = object$sigma, lambda = object$lambda,
                     R2 = object$R2)
  out$weighted <- !is.null(object$weights)
  if (!is.null(object$avgderivatives)) {
    se <- sqrt(as.numeric(object$var.avgderivatives))
    avg <- as.numeric(object$avgderivatives)
    tval <- avg / se
    pval <- 2 * pt(-abs(tval), df = nrow(object$X) - 1L)
    coefs <- cbind(Estimate   = avg,
                   `Std. Err` = se,
                   `t value`  = tval,
                   `Pr(>|t|)` = pval)
    rownames(coefs) <- colnames(object$avgderivatives)
    out$avg_eff <- coefs
    quart <- apply(object$derivatives, 2L, quantile,
                   probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    out$quartiles <- quart
  }
  class(out) <- "summary.krls_rr"
  out
}

#' @export
print.summary.krls_rr <- function(x, ...) {
  cat("Kernel Regularized Least Squares (KRLS)\n")
  ci <- x$call_info
  cat(sprintf("  n = %d  p = %d  sigma = %s  lambda = %s  R^2 = %s\n",
              as.integer(ci["n"]), as.integer(ci["p"]),
              signif(ci["sigma"], 4),
              signif(ci["lambda"], 4),
              signif(ci["R2"], 4)))
  if (isTRUE(x$weighted)) {
    cat("  weighted: yes (mean(w)=1 normalisation)\n")
  }
  if (!is.null(x$avg_eff)) {
    cat("\nAverage Marginal Effects:\n")
    printCoefmat(x$avg_eff, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nQuartiles of Pointwise Marginal Effects:\n")
    print(x$quartiles)
  }
  invisible(x)
}
