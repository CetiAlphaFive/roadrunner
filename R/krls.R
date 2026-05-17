# krls -- Kernel Regularized Least Squares.
#
# File layout:
#   1. krls()              -- main entry; mirrors KRLS::krls (Hainmueller and
#                             Hazlett 2014) numerically.  Returns an S3
#                             object of class "krls".
#   2. predict.krls()      -- predictions on new data, optional fit SEs.
#   3. print.krls()        -- compact one-screen summary.
#   4. summary.krls()      -- richer summary with quartiles of marginal
#                             effects (mirrors KRLS::summary.krls layout).
#   5. .krls_lambdasearch  -- golden-section LOO lambda selector.
#   6. .krls_fd_binary     -- replace binary-column derivatives with
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
#' The numerical pipeline mirrors `KRLS::krls()` exactly:
#'
#' 1. `X` and `y` are standardised (column-centred, unit sd).
#' 2. The Gaussian kernel `K_ij = exp(-||x_i - x_j||^2 / sigma)` is built
#'    in C++ with TBB-parallel inner loops.
#' 3. `K` is eigendecomposed.  All subsequent solves use the
#'    eigen-basis closed forms, never inverting `K + lambda I` directly.
#' 4. `lambda` is selected by golden-section search on the closed-form
#'    LOO error sum, with the same `(L, U, tol)` bracket as `KRLS::krls`.
#' 5. Marginal effects use the identity
#'    `dy/dx_k = -(2/sigma) * (X_k * (K c) - K diag(c) X)_k`, dispatched
#'    through BLAS so no explicit `n x n` distance matrix is materialised.
#' 6. Output is unstandardised back to the original `(X, y)` scale.
#'
#' At a fixed `(sigma, lambda)`, fits agree with `KRLS::krls()` to within
#' floating-point precision (typically `< 1e-10` on coefficients and
#' fitted values for `n <= 1000`).
#'
#' @param X A numeric matrix of predictors (`n x p`).  Constant columns
#'   are rejected.  Missing values are not allowed.
#' @param y A numeric response vector or single-column matrix.  Constant
#'   `y` is rejected.
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
#' @param L,U Optional lower / upper bracket for the lambda search.  If
#'   `NULL`, defaults follow `KRLS::krls()`.
#' @param tol Tolerance for the lambda golden section.  Default
#'   `1e-3 * n`.
#' @param eigtrunc Optional eigenvalue truncation cutoff in `(0, 1]`.
#'   When set, eigenvalues below `eigtrunc * max(d)` are dropped from
#'   the solve.  `NULL` (default) keeps all eigenvalues.
#' @param print.level Integer.  `0` is silent, `1` prints the chosen
#'   lambda and marginal-effect summaries.  Default `0`.
#'
#' @return An object of S3 class `"krls"` with components mirroring
#'   `KRLS::krls()`: `K`, `coeffs`, `Looe`, `fitted`, `X`, `y`, `sigma`,
#'   `lambda`, `R2`, `derivatives`, `avgderivatives`,
#'   `var.avgderivatives`, `vcov.c`, `vcov.fitted`, `binaryindicator`.
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
#' n <- 80
#' x <- matrix(rnorm(n * 2), n, 2)
#' y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.1)
#' fit <- krls(x, y, derivative = TRUE)
#' fit$avgderivatives
#'
#' @export
krls <- function(X, y,
                 sigma = NULL, lambda = NULL,
                 derivative = TRUE, binary = TRUE, vcov = TRUE,
                 L = NULL, U = NULL, tol = NULL, eigtrunc = NULL,
                 print.level = 0) {
  ## --- argument validation -----------------------------------------
  if (is.null(X) || is.null(y)) stop("X and y are required")
  ## Reject non-numeric y BEFORE coercion to avoid silent string->double
  ## via storage.mode (audit KRLS-AUDIT-002).
  if (is.factor(y) || is.character(y)) {
    stop("y must be numeric (got ",
         if (is.factor(y)) "factor" else "character", ")")
  }
  if (is.factor(X) || is.character(X)) {
    stop("X must be numeric (got ",
         if (is.factor(X)) "factor" else "character", ")")
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  ## Reject multi-column y explicitly (audit KRLS-AUDIT-001).  Without
  ## this check the downstream scale() call errors with an opaque
  ## "length of 'center' must equal the number of columns of 'x'".
  if (ncol(y) != 1L) {
    stop("y must be a vector or single-column matrix (got ", ncol(y),
         " columns)")
  }
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  if (!is.numeric(X)) stop("X must be numeric")
  if (!is.numeric(y)) stop("y must be numeric")
  if (anyNA(X))      stop("X contains missing data")
  if (anyNA(y))      stop("y contains missing data")
  if (var(as.vector(y)) == 0) stop("y is a constant (does not vary)")
  stopifnot(is.logical(derivative), is.logical(vcov), is.logical(binary))
  if (derivative && !vcov) {
    stop("derivative = TRUE requires vcov = TRUE")
  }
  n <- nrow(X); d <- ncol(X)
  if (n != nrow(y)) stop("nrow(X) not equal to length of y")
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

  ## --- standardise --------------------------------------------------
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

  ## --- kernel + eigendecomposition ---------------------------------
  K <- krls_kernel_cpp(Xs, sigma)
  eo <- krls_eig_cpp(K)
  dvals <- as.numeric(eo$values)
  V     <- eo$vectors
  Vsq   <- krls_vsq_cpp(V)
  Vty   <- as.numeric(crossprod(V, ys))

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
  if (is.null(lambda)) {
    lambda <- .krls_lambdasearch(dvals = dvals, V = V, Vsq = Vsq,
                                 Vty = Vty, n_y = nrow(y),
                                 L = L, U = U, tol = tol,
                                 noisy = isTRUE(print.level > 2))
    if (print.level > 1) {
      cat("Lambda that minimizes Loo-Loss is:", round(lambda, 5), "\n")
    }
  } else {
    stopifnot(is.numeric(lambda), length(lambda) == 1, lambda > 0)
  }

  ## --- solve for c, LOOe, diag(G^{-1}) -----------------------------
  sol     <- krls_solve_cpp(dvals, V, Vsq, Vty, lambda)
  coeffs  <- as.numeric(sol$coeffs)
  Le      <- sol$Le
  ## yhat on standardised scale
  yfit_s  <- as.numeric(K %*% coeffs)

  ## --- coefficient covariance --------------------------------------
  vcov.c <- NULL
  vcov.fitted <- NULL
  sigma2 <- as.numeric((1 / n) * crossprod(ys - yfit_s))
  if (isTRUE(vcov)) {
    w <- sigma2 / (dvals + lambda)^2
    ## vcov_c_std = V diag(w) V'
    vcovmatc <- V %*% (t(V) * w)
    vcovmatyhat <- crossprod(K, vcovmatc %*% K)
  }

  ## --- derivatives --------------------------------------------------
  derivmat <- avgderiv <- varavgderivmat <- NULL
  if (isTRUE(derivative)) {
    derivmat_s <- krls_deriv_cpp(Xs, K, coeffs, sigma)
    avgderiv_s <- matrix(colMeans(derivmat_s), nrow = 1)
    varavg_s   <- krls_avg_deriv_var_cpp(Xs, K, V, dvals, sigma, lambda, sigma2)
    colnames(derivmat_s)  <- colnames(X)
    colnames(avgderiv_s)  <- colnames(X)
    ## unstandardise: dy/dx_k = (sd_y / sd_xk) * dy_s/dx_s_k
    scale_vec <- as.numeric(y.init.sd) / X.init.sd
    derivmat  <- sweep(derivmat_s, 2L, scale_vec, "*")
    avgderiv  <- matrix(colMeans(derivmat), nrow = 1,
                        dimnames = list(NULL, colnames(X)))
    varavgderivmat <- matrix(scale_vec^2 * as.numeric(varavg_s),
                             nrow = 1,
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
  class(out) <- c("krls_rr", "krls")

  if (isTRUE(derivative) && isTRUE(binary)) {
    out <- .krls_fd_binary(out)
  }

  if (print.level > 0 && isTRUE(derivative)) {
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
#' @param newdata A numeric matrix with the same columns as the training
#'   `X`.
#' @param se.fit Logical.  If `TRUE`, return pointwise standard errors
#'   of the predictions.  Requires the fit was created with
#'   `vcov = TRUE`.
#' @param ... Currently unused.
#' @export
predict.krls_rr <- function(object, newdata, se.fit = FALSE, ...) {
  if (!inherits(object, "krls")) {
    stop("object is not of class 'krls'")
  }
  if (isTRUE(se.fit) && is.null(object$vcov.c)) {
    stop("refit with krls(..., vcov = TRUE) to compute standard errors")
  }
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  ## Reorder newdata columns by name when both training X and newdata
  ## carry colnames (audit KRLS-AUDIT-004).  Otherwise fall back to
  ## positional; warn if the caller passed names and they differ
  ## (mismatch could otherwise be silently wrong).
  cn_train <- colnames(object$X)
  cn_new   <- colnames(newdata)
  if (!is.null(cn_train) && !is.null(cn_new)) {
    if (!setequal(cn_train, cn_new)) {
      stop("colnames(newdata) do not match colnames(X) from fit: ",
           "training = c(", paste(shQuote(cn_train), collapse = ", "),
           "), newdata = c(", paste(shQuote(cn_new), collapse = ", "),
           ")")
    }
    if (!identical(cn_train, cn_new)) {
      newdata <- newdata[, cn_train, drop = FALSE]
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
  Xn     <- scale(newdata, center = Xmeans, scale = Xsd)
  Xn     <- matrix(Xn, nrow(newdata), ncol(newdata))
  nn     <- nrow(Xn)

  Knew   <- krls_kernel_pred_cpp(Xn, Xs, object$sigma)
  yhat   <- Knew %*% object$coeffs

  vcov.fit <- se.fit.out <- NULL
  if (isTRUE(se.fit)) {
    ## vcov.c (stored on raw-y scale) -> standardised-y scale
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
  if (!is.null(x$avg_eff)) {
    cat("\nAverage Marginal Effects:\n")
    printCoefmat(x$avg_eff, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nQuartiles of Pointwise Marginal Effects:\n")
    print(x$quartiles)
  }
  invisible(x)
}
