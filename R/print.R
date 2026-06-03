# S3 print / summary / plot methods for the `ares` class.

#' Print method for `ares` fits
#'
#' Prints the call, family, term count, training RSS and GCV, and the
#' control parameters used. For `family = "binomial"`, also prints the
#' deviance, McFadden pseudo-R^2, AIC, and IRLS convergence.
#'
#' @param x An object of class `"ares"`.
#' @param digits Significant digits for numeric output.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' print(fit)
#' @export
print.ares <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  fam <- if (is.null(x$family)) "gaussian" else x$family
  cat("\nMARS fit:  family =", fam, "\n")
  cat("  Selected terms:", length(x$selected.terms), "of", nrow(x$dirs),
      "forward-pass terms\n")
  cat("  RSS:", format(x$rss, digits = digits), "\n")
  cat("  GCV:", format(x$gcv, digits = digits), "\n")
  if (fam == "binomial" && !is.null(x$glm)) {
    g <- x$glm
    # Pseudo-R^2 (McFadden): 1 - dev / null.dev.
    pr2 <- if (is.finite(g$null.deviance) && g$null.deviance > 0) {
      1 - g$deviance / g$null.deviance
    } else NA_real_
    cat("  Deviance:", format(g$deviance, digits = digits),
        "  Null:", format(g$null.deviance, digits = digits),
        "  McFadden R2:", format(pr2, digits = digits), "\n")
    cat("  AIC:", format(g$aic, digits = digits),
        "  df.resid:", g$df.residual,
        "  IRLS iter:", g$iter,
        "  converged:", g$converged, "\n")
  }
  cat("  Degree:", x$degree, "  Penalty:", x$penalty, "  nthreads:", x$nthreads, "\n")
  # BUG-013 (v0.0.0.9032): surface bagging + autotune state so users can
  # see at a glance that predictions are bag-averaged or that hyper-
  # parameters were cross-validated. Both objects are first-class API
  # surfaces but used to be invisible in print(); the only signal was
  # whatever happened to land in the `Call:` line.
  if (!is.null(x$boot)) {
    cat("  Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  }
  if (!is.null(x$autotune)) {
    at <- x$autotune
    cat(sprintf(
      "  Autotune: degree=%s  penalty=%s  nk=%s  fast.k=%s  warmstart=%s\n",
      format(at$degree), format(at$penalty, digits = digits),
      format(at$nk), format(at$fast_k), format(at$warmstart)))
  }
  invisible(x)
}

#' Summary method for `ares` fits
#'
#' Returns a compact summary of a fitted `ares` model: the call, the
#' table of selected terms and coefficients, training RSS and GCV, and
#' (for non-gaussian families) GLM-side fit statistics.
#'
#' @param object An object of class `"ares"`.
#' @param ... Currently ignored.
#' @return An object of class `"summary.ares"`.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' summary(fit)
#' @export
summary.ares <- function(object, ...) {
  tab <- data.frame(
    term  = names(object$coefficients),
    coef  = unname(object$coefficients),
    stringsAsFactors = FALSE
  )
  fam <- if (is.null(object$family)) "gaussian" else object$family
  out <- list(call = object$call, terms = tab,
              rss = object$rss, gcv = object$gcv,
              n_terms = length(object$selected.terms),
              n_forward = nrow(object$dirs),
              degree = object$degree, penalty = object$penalty,
              family = fam,
              glm = object$glm,
              # BUG-013 (v0.0.0.9032): carry bag + autotune state through
              # to the summary printer.
              boot = object$boot,
              autotune = object$autotune)
  class(out) <- "summary.ares"
  out
}

#' @rdname summary.ares
#' @param x A `summary.ares` object.
#' @export
print.summary.ares <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nMARS summary:  family =", x$family, "\n")
  cat("  Selected:", x$n_terms, "of", x$n_forward, "forward-pass terms\n")
  cat("  RSS =", format(x$rss), "  GCV =", format(x$gcv),
      "  degree =", x$degree, "\n")
  if (identical(x$family, "binomial") && !is.null(x$glm)) {
    g <- x$glm
    pr2 <- if (is.finite(g$null.deviance) && g$null.deviance > 0) {
      1 - g$deviance / g$null.deviance
    } else NA_real_
    cat("  Deviance =", format(g$deviance),
        "  Null =", format(g$null.deviance),
        "  McFadden R2 =", format(pr2), "\n")
    cat("  AIC =", format(g$aic),
        "  df.resid =", g$df.residual,
        "  converged =", g$converged, "\n")
  }
  if (!is.null(x$boot)) {
    cat("  Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  }
  if (!is.null(x$autotune)) {
    at <- x$autotune
    cat(sprintf(
      "  Autotune: degree=%s  penalty=%s  nk=%s  fast.k=%s  warmstart=%s\n",
      format(at$degree), format(at$penalty),
      format(at$nk), format(at$fast_k), format(at$warmstart)))
  }
  cat("\n")
  print(x$terms, row.names = FALSE)
  invisible(x)
}

#' Diagnostic plots for an `ares` fit
#'
#' Six diagnostic panels modelled on [stats::plot.lm()]: residuals vs
#' fitted, normal Q-Q of standardized residuals, scale-location,
#' Cook's distance, residuals vs leverage (with Cook's-distance
#' contours), and Cook's distance vs leverage. The default `which`
#' draws panels 1, 2, 3, 5 in a 2x2 grid when the device is on a
#' single-panel layout; this matches the `plot.lm()` default. For
#' non-gaussian families (`binomial`, `poisson`, `gamma`), deviance
#' residuals are used in panel 1 and IRLS working weights are used to
#' compute leverages and standardized residuals (the GLM analogue of
#' the gaussian hat matrix).
#'
#' @param x A fitted object of class `"ares"`.
#' @param which Integer subset of `1:6` selecting which panels to
#'   draw. Default `c(1, 2, 3, 5)` (the `plot.lm()` default).
#' @param caption Length-6 list of panel captions; defaults match
#'   `plot.lm()`.
#' @param panel Panel function for scatter panels; defaults to
#'   [graphics::panel.smooth()].
#' @param sub.caption Optional caption shown on the outer margin (one
#'   per page). Defaults to a deparsed `x$call`.
#' @param main Title for individual panels (typically `""`).
#' @param ask If `TRUE`, prompt before each page when multiple panels
#'   are drawn on a single-panel device.
#' @param ... Further graphical parameters passed through to the
#'   underlying `plot()` calls.
#' @param id.n Number of extreme points to label per panel (default
#'   `3`; set `0` to suppress).
#' @param labels.id Character vector of row labels; defaults to
#'   `1:n`.
#' @param cex.id `cex` for point labels.
#' @param qqline Logical; if `TRUE`, add a Q-Q reference line on
#'   panel 2.
#' @param cook.levels Cook's-distance contours drawn on panel 5.
#' @param add.smooth Logical; passed to `panel`.
#' @param label.pos `pos` argument forwarded to [graphics::text()].
#' @return Invisibly returns `x`.
#' @examples
#' \dontrun{
#'   fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#'   plot(fit)
#' }
#' @export
plot.ares <- function(x,
                      which = c(1L, 2L, 3L, 5L),
                      caption = list("Residuals vs Fitted", "Normal Q-Q",
                                     "Scale-Location", "Cook's distance",
                                     "Residuals vs Leverage",
                                     expression("Cook's dist vs Leverage  " *
                                                h[ii] / (1 - h[ii]))),
                      panel = if (add.smooth) graphics::panel.smooth
                              else graphics::points,
                      sub.caption = NULL, main = "",
                      ask = prod(graphics::par("mfcol")) <
                              length(which) &&
                            grDevices::dev.interactive(),
                      ...,
                      id.n = 3L,
                      labels.id = NULL,
                      cex.id = 0.75,
                      qqline = TRUE,
                      cook.levels = c(0.5, 1.0),
                      add.smooth = getOption("add.smooth", TRUE),
                      label.pos = c(4, 2)) {
  if (!inherits(x, "ares"))
    stop("plot.ares: 'x' must be an 'ares' object.")
  which <- as.integer(which)
  if (any(!which %in% 1:6))
    stop("plot.ares: `which` must be a subset of 1:6.")
  show <- rep(FALSE, 6L); show[which] <- TRUE

  fam <- if (is.null(x$family)) "gaussian" else x$family
  yhat <- as.numeric(x$fitted.values)
  rraw <- as.numeric(x$residuals)
  n <- length(yhat)
  if (n == 0L) stop("plot.ares: fit has zero observations.")
  w <- if (!is.null(x$weights)) as.numeric(x$weights) else rep(1, n)

  bx <- x$bx
  if (is.null(bx))
    stop("plot.ares: $bx not on fit; diagnostics need the design matrix.")

  # Build sqrt(working-weights) for the hat matrix. Gaussian uses the
  # supplied observation weights; non-gaussian uses canonical-link IRLS
  # weights: W = w * (dmu/deta)^2 / V(mu). On the canonical links
  # implemented here those collapse to mu*(1-mu) (binomial),
  # mu (poisson), 1 (gamma), times the supplied `w`.
  if (fam == "gaussian") {
    sqW <- sqrt(w)
    r_for_plot <- rraw
    # Recover the response so we can also compute a Pearson-style
    # standardized resid below in the same code path. y = yhat + r_raw.
    yv <- yhat + rraw
  } else {
    mu <- yhat
    yv <- yhat + rraw
    Wfac <- switch(fam,
                   binomial = pmax(mu * (1 - mu), 1e-12),
                   poisson  = pmax(mu, 1e-12),
                   gamma    = rep(1, n))
    W <- w * Wfac
    sqW <- sqrt(pmax(W, 0))
    fam_obj <- switch(fam,
                      binomial = stats::binomial(),
                      poisson  = stats::poisson(),
                      gamma    = stats::Gamma())
    d_i <- fam_obj$dev.resids(yv, mu, w)
    r_for_plot <- sign(yv - mu) * sqrt(pmax(d_i, 0))
  }

  # Hat values via QR of sqrt(W) * bx. rkB is the effective rank used
  # in standardized-resid and Cook's denominators.
  Btil <- sqW * bx
  qrB <- qr(Btil)
  rkB <- max(qrB$rank, 1L)
  Q <- qr.Q(qrB)
  if (ncol(Q) > rkB) Q <- Q[, seq_len(rkB), drop = FALSE]
  h <- rowSums(Q * Q)
  h <- pmin(pmax(h, 0), 1 - .Machine$double.eps)

  # Dispersion + standardized residuals. plot.lm uses
  # rstandard = r * sqrt(w) / (sigma * sqrt(1 - h)). For binomial /
  # poisson the dispersion is fixed at 1; for gamma we estimate phi
  # from the Pearson residuals on the selected basis.
  if (fam == "gaussian") {
    rdf <- max(n - rkB, 1L)
    sigma <- sqrt(sum(w * rraw^2) / rdf)
    rstd <- rraw * sqrt(w) / (sigma * sqrt(1 - h))
  } else if (fam == "gamma") {
    pear <- (yv - yhat) / sqrt(yhat^2 / pmax(w, 1e-12))
    rdf <- max(n - rkB, 1L)
    phi <- sum(pear^2) / rdf
    rstd <- r_for_plot / sqrt(phi * (1 - h))
  } else {
    rstd <- r_for_plot / sqrt(1 - h)
  }
  rstd[!is.finite(rstd)] <- NA_real_

  cook <- (rstd^2 / rkB) * h / (1 - h)
  cook[!is.finite(cook)] <- NA_real_

  if (is.null(labels.id)) labels.id <- as.character(seq_len(n))
  extrm <- function(v, k = id.n) {
    if (k < 1L) return(integer(0))
    finite <- which(is.finite(v))
    if (!length(finite)) return(integer(0))
    ord <- finite[order(-abs(v[finite]))]
    ord[seq_len(min(k, length(ord)))]
  }

  # Auto-2x2 layout for the natural "dump the diagnostics" call. Only
  # activates when the device is on the default single-panel layout
  # AND multiple panels were requested. plot.lm leaves this to the
  # caller; here we default to the convenience layout.
  one_fig <- all(graphics::par("mfcol") == c(1L, 1L))
  if (one_fig && length(which) > 1L) {
    op <- graphics::par(mfrow = c(2L, 2L), oma = c(0, 0, 2, 0),
                        mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
  } else if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op), add = TRUE)
  }

  if (is.null(sub.caption)) {
    cl <- if (!is.null(x$call))
      paste(deparse(x$call, width.cutoff = 75L), collapse = " ")
    else "ares fit"
    sub.caption <- if (nchar(cl) > 90L)
      paste0(substring(cl, 1L, 87L), "...") else cl
  }
  resid_lbl <- if (fam == "gaussian") "Residuals" else "Deviance residuals"
  fit_lbl   <- if (fam == "gaussian") "Fitted values"
               else paste0("Predicted values (", fam, ")")

  # Panel 1: residuals vs fitted
  if (show[1L]) {
    ylim <- range(r_for_plot, na.rm = TRUE)
    if (id.n > 0L) ylim <- grDevices::extendrange(r = ylim, f = 0.08)
    graphics::plot(yhat, r_for_plot, xlab = fit_lbl, ylab = resid_lbl,
                   main = main, ylim = ylim, type = "n", ...)
    panel(yhat, r_for_plot, ...)
    graphics::abline(h = 0, lty = 3, col = "gray")
    idx <- extrm(r_for_plot)
    if (length(idx))
      graphics::text(yhat[idx], r_for_plot[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[1L]], side = 3, line = 0.25, cex = 0.85)
  }
  # Panel 2: normal Q-Q of standardized residuals
  if (show[2L]) {
    qq <- stats::qqnorm(rstd, main = main,
                        ylab = "Standardized residuals", ...)
    if (qqline) stats::qqline(rstd, lty = 3, col = "gray")
    idx <- extrm(rstd)
    if (length(idx))
      graphics::text(qq$x[idx], qq$y[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[2L]], side = 3, line = 0.25, cex = 0.85)
  }
  # Panel 3: scale-location
  if (show[3L]) {
    sqrtabs <- sqrt(abs(rstd))
    graphics::plot(yhat, sqrtabs, xlab = fit_lbl,
                   ylab = expression(sqrt(abs(`Standardized residuals`))),
                   main = main, type = "n", ...)
    panel(yhat, sqrtabs, ...)
    idx <- extrm(sqrtabs)
    if (length(idx))
      graphics::text(yhat[idx], sqrtabs[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[3L]], side = 3, line = 0.25, cex = 0.85)
  }
  # Panel 4: Cook's distance
  if (show[4L]) {
    ylim <- c(0, max(cook, na.rm = TRUE) * 1.075)
    graphics::plot(seq_len(n), cook, type = "h", main = main,
                   xlab = "Obs. number", ylab = "Cook's distance",
                   ylim = ylim, ...)
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(idx, cook[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[4L]], side = 3, line = 0.25, cex = 0.85)
  }
  # Panel 5: standardized residuals vs leverage with Cook contours
  if (show[5L]) {
    xlim <- c(0, max(h, na.rm = TRUE) * 1.05)
    ylim <- range(rstd, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-3, 3)
    graphics::plot(h, rstd, xlim = xlim, ylim = ylim, xlab = "Leverage",
                   ylab = "Standardized residuals", main = main,
                   type = "n", ...)
    panel(h, rstd, ...)
    graphics::abline(h = 0, v = 0, lty = 3, col = "gray")
    if (length(cook.levels)) {
      hh <- seq.int(0.001, max(h, na.rm = TRUE), length.out = 101L)
      hh <- hh[hh < 1 & hh > 0]
      for (cl_lvl in cook.levels) {
        rr <- sqrt(cl_lvl * rkB * (1 - hh) / hh)
        graphics::lines(hh, rr, lty = 2, col = "red")
        graphics::lines(hh, -rr, lty = 2, col = "red")
      }
    }
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(h[idx], rstd[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[5L]], side = 3, line = 0.25, cex = 0.85)
  }
  # Panel 6: Cook's vs h/(1-h)
  if (show[6L]) {
    hsc <- h / (1 - h)
    graphics::plot(hsc, cook,
                   xlab = expression(h[ii] / (1 - h[ii])),
                   ylab = "Cook's distance", main = main, type = "n", ...)
    panel(hsc, cook, ...)
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(hsc[idx], cook[idx], labels = labels.id[idx],
                     cex = cex.id, pos = label.pos[1L])
    graphics::mtext(caption[[6L]], side = 3, line = 0.25, cex = 0.85)
  }

  if (!is.null(sub.caption) && length(which) > 1L && one_fig)
    graphics::mtext(sub.caption, outer = TRUE, cex = 0.9, line = 0.5)

  invisible(x)
}

#' Print method for `plda` fits
#'
#' @param x A `"plda"` object.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @examples
#' fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
#' print(fit)
#' @export
print.plda <- function(x, ...) {
  cat("Penalized LDA (", x$penalty, " penalty)\n", sep = "")
  cat("  Classes:      ", paste(x$classes, collapse = ", "), "\n", sep = "")
  cat("  Discriminants:", x$K, "\n")
  cat("  lambda:       ", format(x$lambda, digits = 4), "\n", sep = "")
  if (x$penalty == "fused")
    cat("  lambda2:      ", format(x$lambda2, digits = 4), "\n", sep = "")
  nz <- colSums(abs(x$discrim) > 0)
  cat("  Nonzero feats per discriminant:", paste(nz, collapse = ", "), "\n")
  invisible(x)
}

#' Summary method for `plda` fits
#'
#' @param object A `"plda"` object.
#' @param ... Currently ignored.
#' @return An object of class `"summary.plda"`.
#' @examples
#' fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
#' summary(fit)
#' @export
summary.plda <- function(object, ...) {
  nz <- colSums(abs(object$discrim) > 0)
  structure(list(penalty  = object$penalty,
                 classes  = object$classes,
                 K        = object$K,
                 lambda   = object$lambda,
                 lambda2  = object$lambda2,
                 nonzero  = nz,
                 npred    = nrow(object$discrim),
                 cv       = object$cv),
            class = "summary.plda")
}

#' @rdname summary.plda
#' @param x A `summary.plda` object.
#' @return Invisibly returns `x`.
#' @examples
#' fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
#' print(summary(fit))
#' @export
print.summary.plda <- function(x, ...) {
  cat("Penalized LDA summary\n")
  cat("  penalty:", x$penalty, "| classes:", length(x$classes),
      "| discriminants:", x$K, "\n")
  cat("  lambda:", format(x$lambda, digits = 4), "\n")
  for (k in seq_len(x$K))
    cat(sprintf("  discriminant %d: %d / %d nonzero features\n",
                k, x$nonzero[k], x$npred))
  if (!is.null(x$cv))
    cat("  CV-selected via", length(x$cv$grid), "lambda grid points\n")
  invisible(x)
}

#' Projection plot for `plda` fits
#'
#' @param x A `"plda"` object.
#' @param data Predictor matrix (required).
#' @param labels Optional class label vector for colouring points.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return Invisibly returns `x`.
#' @examples
#' \dontrun{
#'   fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
#'   plot(fit, data = iris[, 1:4])
#' }
#' @export
plot.plda <- function(x, data = NULL, labels = NULL, ...) {
  if (is.null(data))
    stop("plot.plda: supply `data` (predictor matrix used to fit).",
         call. = FALSE)
  sc <- plda_project_cpp(as.matrix(data), x$mu, x$sdw, x$discrim)
  cols <- if (is.null(labels)) 1L else as.integer(as.factor(labels))
  if (x$K >= 2L) {
    graphics::plot(sc[, 1], sc[, 2], col = cols, pch = 19,
                   xlab = "Discriminant 1", ylab = "Discriminant 2",
                   main = "Penalized LDA projection", ...)
  } else {
    d <- sc[, 1]
    graphics::plot(d, jitter(rep(0, length(d))),
                   col = cols, pch = 19, yaxt = "n",
                   xlab = "Discriminant 1", ylab = "",
                   main = "Penalized LDA projection", ...)
  }
  invisible(x)
}

# ============================================================================
#  S3 print / summary / plot methods for the `ols` class
# ============================================================================

#' Print method for `ols` fits
#'
#' Prints the call, the fitted coefficients, the residual standard error,
#' and the residual degrees of freedom.
#'
#' @param x An object of class `"ols"`.
#' @param digits Significant digits for numeric output.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @examples
#' fit <- ols(mpg ~ wt + hp, data = mtcars)
#' print(fit)
#' @export
print.ols <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  weighted <- !is.null(x$weights)
  cat("\n", if (weighted) "Weighted " else "",
      "linear model (OLS)\n", sep = "")
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\n  Residual SE:", format(x$sigma, digits = digits),
      "on", x$df.residual, "degrees of freedom\n")
  if (!is.null(x$boot))
    cat("  Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  if (!is.null(x$varmod))
    cat("  Variance model:", x$varmod$type, "\n")
  invisible(x)
}

#' Summary method for `ols` fits
#'
#' Returns the coefficient table (estimate, standard error, t value,
#' p value), the residual standard error, R-squared and adjusted
#' R-squared, and the overall F-statistic.
#'
#' @param object An object of class `"ols"`.
#' @param robust Heteroscedasticity-consistent covariance for the
#'   standard errors and t / F tests. `"none"` (default) uses the
#'   classical covariance; `"HC0"`-`"HC3"` use a sandwich estimator.
#' @param ... Currently ignored.
#' @return An object of class `"summary.ols"`.
#' @examples
#' fit <- ols(mpg ~ wt + hp, data = mtcars)
#' summary(fit)
#' summary(fit, robust = "HC3")
#' @export
summary.ols <- function(object,
                        robust = c("none", "HC0", "HC1", "HC2", "HC3"),
                        ...) {
  robust <- match.arg(robust)
  vcov <- .ols_vcov(object, robust)
  est <- object$coefficients
  se <- sqrt(diag(vcov))
  tval <- est / se
  df <- object$df.residual
  pval <- 2 * stats::pt(-abs(tval), df = df)
  coef_tab <- cbind(Estimate = est, `Std. Error` = se,
                    `t value` = tval, `Pr(>|t|)` = pval)
  rownames(coef_tab) <- names(est)

  # R-squared. With weights, the total sum of squares is weighted around
  # the weighted mean; the intercept governs whether the mean is removed.
  w <- if (is.null(object$weights)) rep(1, object$n) else
    as.numeric(object$weights)
  y <- object$y
  has_int <- isTRUE(object$intercept) ||
    "(Intercept)" %in% names(est)
  ybar <- if (has_int) stats::weighted.mean(y, w) else 0
  tss <- sum(w * (y - ybar)^2)
  rss <- object$rss
  r_squared <- if (tss > 0) 1 - rss / tss else NA_real_
  # df for the model: p minus the intercept when present.
  p_model <- object$p - as.integer(has_int)
  adj_r_squared <- if (tss > 0 && df > 0)
    1 - (1 - r_squared) * (object$n - as.integer(has_int)) / df else
      NA_real_

  # Overall F-statistic: model mean square over residual mean square.
  # When robust != "none" the test is a Wald F on the non-intercept
  # coefficients using the robust covariance; otherwise the classical F.
  fstat <- NULL
  if (p_model >= 1L) {
    if (identical(robust, "none")) {
      f_value <- ((tss - rss) / p_model) / (rss / df)
    } else {
      idx <- if (has_int) which(names(est) != "(Intercept)") else
        seq_along(est)
      b <- est[idx]
      vb <- vcov[idx, idx, drop = FALSE]
      wald <- as.numeric(t(b) %*% solve(vb, b))
      f_value <- wald / p_model
    }
    fstat <- list(value = f_value, numdf = p_model, dendf = df,
                  pvalue = stats::pf(f_value, p_model, df,
                                     lower.tail = FALSE))
  }

  structure(list(call = object$call, coefficients = coef_tab,
                 sigma = object$sigma, df.residual = df,
                 r.squared = r_squared, adj.r.squared = adj_r_squared,
                 fstatistic = fstat, robust = robust,
                 weighted = !is.null(object$weights),
                 boot = object$boot),
            class = "summary.ols")
}

#' @rdname summary.ols
#' @param x A `summary.ols` object.
#' @param digits Significant digits for numeric output.
#' @export
print.summary.ols <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("Call:\n"); print(x$call)
  cat("\n", if (x$weighted) "Weighted " else "",
      "linear model (OLS)", sep = "")
  if (!identical(x$robust, "none"))
    cat("  --  ", x$robust, " robust standard errors", sep = "")
  cat("\n\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = getOption("show.signif.stars", TRUE),
                      has.Pvalue = TRUE)
  cat("\nResidual standard error:", format(x$sigma, digits = digits),
      "on", x$df.residual, "degrees of freedom\n")
  cat("Multiple R-squared: ", format(x$r.squared, digits = digits),
      ",  Adjusted R-squared: ", format(x$adj.r.squared, digits = digits),
      "\n", sep = "")
  if (!is.null(x$fstatistic)) {
    f <- x$fstatistic
    cat("F-statistic: ", format(f$value, digits = digits),
        " on ", f$numdf, " and ", f$dendf, " DF,  p-value: ",
        format.pval(f$pvalue, digits = digits), "\n", sep = "")
  }
  if (!is.null(x$boot))
    cat("Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  invisible(x)
}

#' Diagnostic plots for an `ols` fit
#'
#' Four diagnostic panels modelled on [stats::plot.lm()]: residuals vs
#' fitted, normal Q-Q of standardized residuals, scale-location, and
#' residuals vs leverage with Cook's-distance contours.
#'
#' @param x A fitted object of class `"ols"`.
#' @param which Integer subset of `1:4` selecting which panels to draw.
#'   Default `1:4`.
#' @param id.n Number of extreme points to label per panel (default
#'   `3`; set `0` to suppress).
#' @param ... Further graphical parameters passed to the underlying
#'   `plot()` calls.
#' @return Invisibly returns `x`.
#' @examples
#' \dontrun{
#'   fit <- ols(mpg ~ wt + hp, data = mtcars)
#'   plot(fit)
#' }
#' @export
plot.ols <- function(x, which = 1:4, id.n = 3L, ...) {
  if (!inherits(x, "ols"))
    stop("plot.ols: 'x' must be an 'ols' object.")
  which <- as.integer(which)
  if (any(!which %in% 1:4))
    stop("plot.ols: `which` must be a subset of 1:4.")
  show <- rep(FALSE, 4L); show[which] <- TRUE

  yhat <- as.numeric(x$fitted.values)
  rraw <- as.numeric(x$residuals)
  n <- length(yhat)
  if (n == 0L) stop("plot.ols: fit has zero observations.")
  w <- if (!is.null(x$weights)) as.numeric(x$weights) else rep(1, n)
  h <- pmin(pmax(as.numeric(x$hatvalues), 0), 1 - .Machine$double.eps)
  rk <- max(x$rank, 1L)
  sigma <- x$sigma
  # plot.lm standardized residual: r sqrt(w) / (sigma sqrt(1 - h)).
  rstd <- rraw * sqrt(w) / (sigma * sqrt(1 - h))
  rstd[!is.finite(rstd)] <- NA_real_
  cook <- (rstd^2 / rk) * h / (1 - h)
  cook[!is.finite(cook)] <- NA_real_

  labels.id <- as.character(seq_len(n))
  extrm <- function(v, k = id.n) {
    if (k < 1L) return(integer(0))
    finite <- which(is.finite(v))
    if (!length(finite)) return(integer(0))
    ord <- finite[order(-abs(v[finite]))]
    ord[seq_len(min(k, length(ord)))]
  }

  one_fig <- all(graphics::par("mfcol") == c(1L, 1L))
  if (one_fig && length(which) > 1L) {
    op <- graphics::par(mfrow = c(2L, 2L), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  # Panel 1: residuals vs fitted.
  if (show[1L]) {
    graphics::plot(yhat, rraw, xlab = "Fitted values", ylab = "Residuals",
                   main = "Residuals vs Fitted", ...)
    graphics::abline(h = 0, lty = 3, col = "gray")
    graphics::lines(stats::lowess(yhat, rraw), col = "red")
    idx <- extrm(rraw)
    if (length(idx))
      graphics::text(yhat[idx], rraw[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 2: normal Q-Q of standardized residuals.
  if (show[2L]) {
    qq <- stats::qqnorm(rstd, main = "Normal Q-Q",
                        ylab = "Standardized residuals", ...)
    stats::qqline(rstd, lty = 3, col = "gray")
    idx <- extrm(rstd)
    if (length(idx))
      graphics::text(qq$x[idx], qq$y[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 3: scale-location.
  if (show[3L]) {
    sqrtabs <- sqrt(abs(rstd))
    graphics::plot(yhat, sqrtabs, xlab = "Fitted values",
                   ylab = expression(sqrt(abs(`Standardized residuals`))),
                   main = "Scale-Location", ...)
    ok <- is.finite(yhat) & is.finite(sqrtabs)
    if (any(ok))
      graphics::lines(stats::lowess(yhat[ok], sqrtabs[ok]), col = "red")
    idx <- extrm(sqrtabs)
    if (length(idx))
      graphics::text(yhat[idx], sqrtabs[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 4: residuals vs leverage with Cook's-distance contours.
  if (show[4L]) {
    xlim <- c(0, max(h, na.rm = TRUE) * 1.05)
    ylim <- range(rstd, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-3, 3)
    graphics::plot(h, rstd, xlim = xlim, ylim = ylim,
                   xlab = "Leverage", ylab = "Standardized residuals",
                   main = "Residuals vs Leverage", ...)
    graphics::abline(h = 0, v = 0, lty = 3, col = "gray")
    hh <- seq.int(0.001, max(h, na.rm = TRUE), length.out = 101L)
    hh <- hh[hh < 1 & hh > 0]
    for (cl_lvl in c(0.5, 1.0)) {
      rr <- sqrt(cl_lvl * rk * (1 - hh) / hh)
      graphics::lines(hh, rr, lty = 2, col = "red")
      graphics::lines(hh, -rr, lty = 2, col = "red")
    }
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(h[idx], rstd[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  invisible(x)
}

# ============================================================================
#  S3 print / summary / plot methods for the `logreg` class
# ============================================================================

#' Print method for `logreg` fits
#'
#' Prints the call, the fitted coefficients, the residual and null
#' deviance, the IRLS iteration count, and the convergence status.
#'
#' @param x An object of class `"logreg"`.
#' @param digits Significant digits for numeric output.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @examples
#' df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
#' print(logreg(y ~ wt + hp, data = df))
#' @export
print.logreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  cat("Call:\n")
  print(x$call)
  weighted <- !is.null(x$weights)
  cat("\n", if (weighted) "Weighted " else "",
      "binary logistic regression\n", sep = "")
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\n  Null deviance:    ", format(x$null.deviance, digits = digits),
      " on ", x$df.null, " degrees of freedom\n", sep = "")
  cat("  Residual deviance:", format(x$deviance, digits = digits),
      " on ", x$df.residual, " degrees of freedom\n", sep = "")
  cat("  AIC:", format(x$aic, digits = digits), "\n")
  cat("  IRLS iterations:", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n")
  if (!is.null(x$boot))
    cat("  Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  invisible(x)
}

#' Summary method for `logreg` fits
#'
#' Returns the coefficient table (estimate, standard error, z value,
#' p value), the residual and null deviance, the AIC, and the IRLS
#' iteration count / convergence status.
#'
#' @param object An object of class `"logreg"`.
#' @param robust Heteroscedasticity-consistent covariance for the
#'   standard errors and z tests. `"none"` (default) uses the classical
#'   maximum-likelihood covariance; `"HC0"`-`"HC3"` use a sandwich
#'   estimator.
#' @param ... Currently ignored.
#' @return An object of class `"summary.logreg"`.
#' @examples
#' df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
#' fit <- logreg(y ~ wt + hp, data = df)
#' summary(fit)
#' summary(fit, robust = "HC3")
#' @export
summary.logreg <- function(object,
                           robust = c("none", "HC0", "HC1", "HC2", "HC3"),
                           ...) {
  robust <- match.arg(robust)
  vcov <- .logreg_vcov(object, robust)
  est <- object$coefficients
  se <- sqrt(diag(vcov))
  zval <- est / se
  pval <- 2 * stats::pnorm(-abs(zval))
  coef_tab <- cbind(Estimate = est, `Std. Error` = se,
                    `z value` = zval, `Pr(>|z|)` = pval)
  rownames(coef_tab) <- names(est)

  structure(list(call = object$call, coefficients = coef_tab,
                 deviance = object$deviance,
                 null.deviance = object$null.deviance,
                 df.residual = object$df.residual,
                 df.null = object$df.null,
                 aic = object$aic, iter = object$iter,
                 converged = object$converged, robust = robust,
                 weighted = !is.null(object$weights),
                 boot = object$boot),
            class = "summary.logreg")
}

#' @rdname summary.logreg
#' @param x A `summary.logreg` object.
#' @param digits Significant digits for numeric output.
#' @export
print.summary.logreg <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("Call:\n"); print(x$call)
  cat("\n", if (x$weighted) "Weighted " else "",
      "binary logistic regression", sep = "")
  if (!identical(x$robust, "none"))
    cat("  --  ", x$robust, " robust standard errors", sep = "")
  cat("\n\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = getOption("show.signif.stars", TRUE),
                      has.Pvalue = TRUE)
  cat("\n    Null deviance: ", format(x$null.deviance, digits = digits),
      "  on ", x$df.null, " degrees of freedom\n", sep = "")
  cat("Residual deviance: ", format(x$deviance, digits = digits),
      "  on ", x$df.residual, " degrees of freedom\n", sep = "")
  cat("AIC: ", format(x$aic, digits = digits), "\n", sep = "")
  cat("Number of Fisher Scoring iterations: ", x$iter,
      if (x$converged) "" else "  (NOT converged)", "\n", sep = "")
  if (!is.null(x$boot))
    cat("Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")
  invisible(x)
}

#' Diagnostic plots for a `logreg` fit
#'
#' Four diagnostic panels for a logistic-regression fit: deviance
#' residuals vs the linear predictor, a normal Q-Q plot of the deviance
#' residuals, a scale-location panel, and deviance residuals vs leverage
#' with Cook's-distance contours.
#'
#' @param x A fitted object of class `"logreg"`.
#' @param which Integer subset of `1:4` selecting which panels to draw.
#'   Default `1:4`.
#' @param id.n Number of extreme points to label per panel (default
#'   `3`; set `0` to suppress).
#' @param ... Further graphical parameters passed to the underlying
#'   `plot()` calls.
#' @return Invisibly returns `x`.
#' @examples
#' \dontrun{
#'   df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
#'   plot(logreg(y ~ wt + hp, data = df))
#' }
#' @export
plot.logreg <- function(x, which = 1:4, id.n = 3L, ...) {
  if (!inherits(x, "logreg"))
    stop("plot.logreg: 'x' must be a 'logreg' object.")
  which <- as.integer(which)
  if (any(!which %in% 1:4))
    stop("plot.logreg: `which` must be a subset of 1:4.")
  show <- rep(FALSE, 4L); show[which] <- TRUE

  eta <- as.numeric(x$linear.predictors)
  rdev <- as.numeric(x$residuals)              # deviance residuals.
  n <- length(eta)
  if (n == 0L) stop("plot.logreg: fit has zero observations.")
  h <- pmin(pmax(as.numeric(x$hatvalues), 0), 1 - .Machine$double.eps)
  rk <- max(x$rank, 1L)
  # Standardized deviance residual: r_dev / sqrt(1 - h).
  rstd <- rdev / sqrt(1 - h)
  rstd[!is.finite(rstd)] <- NA_real_
  # Approximate Cook's distance for a GLM.
  pear <- as.numeric((x$y - x$fitted.values) /
                       sqrt(x$fitted.values * (1 - x$fitted.values)))
  cook <- (pear^2 / rk) * h / (1 - h)^2
  cook[!is.finite(cook)] <- NA_real_

  labels.id <- as.character(seq_len(n))
  extrm <- function(v, k = id.n) {
    if (k < 1L) return(integer(0))
    finite <- which(is.finite(v))
    if (!length(finite)) return(integer(0))
    ord <- finite[order(-abs(v[finite]))]
    ord[seq_len(min(k, length(ord)))]
  }

  one_fig <- all(graphics::par("mfcol") == c(1L, 1L))
  if (one_fig && length(which) > 1L) {
    op <- graphics::par(mfrow = c(2L, 2L), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  # Panel 1: deviance residuals vs the linear predictor.
  if (show[1L]) {
    graphics::plot(eta, rdev, xlab = "Linear predictor",
                   ylab = "Deviance residuals",
                   main = "Residuals vs Linear Predictor", ...)
    graphics::abline(h = 0, lty = 3, col = "gray")
    ok <- is.finite(eta) & is.finite(rdev)
    if (sum(ok) > 1L)
      graphics::lines(stats::lowess(eta[ok], rdev[ok]), col = "red")
    idx <- extrm(rdev)
    if (length(idx))
      graphics::text(eta[idx], rdev[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 2: normal Q-Q of standardized deviance residuals.
  if (show[2L]) {
    qq <- stats::qqnorm(rstd, main = "Normal Q-Q",
                        ylab = "Std. deviance residuals", ...)
    stats::qqline(rstd, lty = 3, col = "gray")
    idx <- extrm(rstd)
    if (length(idx))
      graphics::text(qq$x[idx], qq$y[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 3: scale-location.
  if (show[3L]) {
    sqrtabs <- sqrt(abs(rstd))
    graphics::plot(eta, sqrtabs, xlab = "Linear predictor",
                   ylab = expression(sqrt(abs(`Std. deviance residuals`))),
                   main = "Scale-Location", ...)
    ok <- is.finite(eta) & is.finite(sqrtabs)
    if (sum(ok) > 1L)
      graphics::lines(stats::lowess(eta[ok], sqrtabs[ok]), col = "red")
    idx <- extrm(sqrtabs)
    if (length(idx))
      graphics::text(eta[idx], sqrtabs[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  # Panel 4: deviance residuals vs leverage with Cook's-distance contours.
  if (show[4L]) {
    xlim <- c(0, max(h, na.rm = TRUE) * 1.05)
    ylim <- range(rstd, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(-3, 3)
    graphics::plot(h, rstd, xlim = xlim, ylim = ylim,
                   xlab = "Leverage", ylab = "Std. deviance residuals",
                   main = "Residuals vs Leverage", ...)
    graphics::abline(h = 0, v = 0, lty = 3, col = "gray")
    hmax <- max(h, na.rm = TRUE)
    if (is.finite(hmax) && hmax > 0) {
      hh <- seq.int(0.001, hmax, length.out = 101L)
      hh <- hh[hh < 1 & hh > 0]
      for (cl_lvl in c(0.5, 1.0)) {
        rr <- sqrt(cl_lvl * rk * (1 - hh) / hh)
        graphics::lines(hh, rr, lty = 2, col = "red")
        graphics::lines(hh, -rr, lty = 2, col = "red")
      }
    }
    idx <- extrm(cook)
    if (length(idx))
      graphics::text(h[idx], rstd[idx], labels = labels.id[idx],
                     cex = 0.75, pos = 4)
  }
  invisible(x)
}

# ============================================================================
#  S3 print / summary / plot methods for the `bgam` class
# ============================================================================

#' Print a `bgam` fit
#'
#' Prints a compact summary of a fitted component-wise P-spline boosting
#' model: the call, family, sample size, number of base-learners, boosting
#' parameters (`mstop`, `nu`, `nknots`, `degree`, `dpen`, `df_target`),
#' CV-selected `mstop_opt` (when `autotune = TRUE`), bagging info (when
#' `n.boot > 0`), the top-5 predictors ranked by selection frequency, and
#' the residual sigma for gaussian fits.
#'
#' @param x An object of class `"bgam"`, as returned by [bgam()].
#' @param digits Number of significant digits for numeric output.  Default
#'   `max(3L, getOption("digits") - 3L)`.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @seealso [bgam()], [summary.bgam()], [plot.bgam()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 4), n, 4)
#' y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
#' fit <- bgam(x, y, mstop = 50, autotune = FALSE)
#' print(fit)
#' @export
print.bgam <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nbgam fit:  family =", x$family, "\n")
  cat("  n =", x$n, "  p (predictors) =", x$p, "\n")
  cat("  mstop =", x$mstop, "  nu =", format(x$nu, digits = digits), "\n")
  cat("  nknots =", x$nknots, "  degree =", x$degree,
      "  dpen =", x$dpen, "  df_target =", x$df_target, "\n")
  if (!is.null(x$cv))
    cat("  Autotune: mstop_opt =", x$cv$mstop_opt, "\n")
  if (!is.null(x$boot))
    cat("  Bagging: n.boot =", x$boot$n.boot, "replicate(s)\n")

  sf <- sort(x$selection_frequency, decreasing = TRUE)
  top5 <- utils::head(sf, 5L)
  cat("\nTop-5 predictors by selection frequency:\n")
  for (nm in names(top5))
    cat("  ", format(nm, width = 12), ":",
        format(top5[nm], digits = digits), "\n")

  if (identical(x$family, "gaussian") && !is.null(x$sigma2))
    cat("\nResidual sigma:", format(sqrt(x$sigma2), digits = digits), "\n")

  invisible(x)
}

#' Summarise a `bgam` fit
#'
#' Returns a `summary.bgam` S3 object containing the call, family,
#' sample and predictor dimensions, boosting parameters (`mstop`, `nu`,
#' `df_target`), CV-selected `mstop_opt` (when `autotune = TRUE`), bagging
#' info (when `n.boot > 0`), training loss at `mstop`, residual sigma for
#' gaussian fits, and a full predictor selection-frequency table sorted
#' in descending order.  The selection frequency is the fraction of
#' boosting iterations in which each predictor was the active base-learner.
#'
#' @param object An object of class `"bgam"`, as returned by [bgam()].
#' @param ... Currently ignored.
#' @return An object of class `"summary.bgam"`.  Print it with
#'   [print.summary.bgam()].
#' @seealso [bgam()], [print.bgam()], [plot.bgam()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 4), n, 4)
#' y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
#' fit <- bgam(x, y, mstop = 50, autotune = FALSE)
#' summary(fit)
#' @export
summary.bgam <- function(object, ...) {
  sf <- sort(object$selection_frequency, decreasing = TRUE)
  sel_tab <- data.frame(
    predictor      = names(sf),
    sel_frequency  = unname(sf),
    stringsAsFactors = FALSE
  )
  out <- list(
    call          = object$call,
    family        = object$family,
    n             = object$n,
    p             = object$p,
    mstop         = object$mstop,
    nu            = object$nu,
    df_target     = object$df_target,
    sel_tab       = sel_tab,
    loss_at_mstop = utils::tail(object$loss_path, 1L),
    sigma2        = object$sigma2,
    cv            = object$cv,
    boot          = object$boot
  )
  class(out) <- "summary.bgam"
  out
}

#' @rdname summary.bgam
#' @param x A `summary.bgam` object returned by [summary.bgam()].
#' @param digits Number of significant digits for numeric output.  Default
#'   `max(3L, getOption("digits") - 3L)`.
#' @return Invisibly returns `x`.
#' @export
print.summary.bgam <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("Call:\n"); print(x$call)
  cat("\nbgam summary:  family =", x$family, "\n")
  cat("  n =", x$n, "  p =", x$p, "\n")
  cat("  mstop =", x$mstop, "  nu =", format(x$nu, digits = digits),
      "  df_target =", x$df_target, "\n")
  if (!is.null(x$cv))
    cat("  CV mstop_opt =", x$cv$mstop_opt, "\n")
  if (!is.null(x$boot))
    cat("  Bagging: n.boot =", x$boot$n.boot, "\n")
  metric <- if (identical(x$family, "binomial")) "deviance" else "mse"
  cat("  Training", metric, "at mstop:",
      format(x$loss_at_mstop, digits = digits), "\n")
  if (!is.null(x$sigma2))
    cat("  Residual sigma:", format(sqrt(x$sigma2), digits = digits), "\n")

  cat("\nPredictor selection frequencies (sorted):\n")
  print(x$sel_tab, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Plot a `bgam` fit
#'
#' Produces a 4-panel diagnostic display for a fitted component-wise P-spline
#' boosting model:
#'
#' \enumerate{
#'   \item Partial effect curve \eqn{\hat{f}_{j^*}(x_{j^*})}{f-hat_j*(x_j*)}
#'     for the top predictor by selection frequency, plotted over a fine grid
#'     covering the training range.
#'   \item Selection-frequency bar chart for all predictors (capped at 20
#'     bars), sorted descending.  Use this panel to identify which predictors
#'     drive the model.
#'   \item Training loss path (MSE for gaussian, mean deviance per observation
#'     for binomial) vs boosting iteration.  A red dashed vertical line marks
#'     `mstop_opt` when `autotune = TRUE`.
#'   \item Fitted vs Observed (gaussian) or fitted-probability box plots by
#'     class (binomial).
#' }
#'
#' The function works in headless (non-interactive) mode without error.
#'
#' @param x An object of class `"bgam"`, as returned by [bgam()].
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @seealso [bgam()], [print.bgam()], [summary.bgam()]
#' @references
#' Buehlmann, P. and Yu, B. (2003). Boosting with the L2 loss: Regression and
#' classification. *Journal of the American Statistical Association*,
#' 98(462):324--339. \doi{10.1198/016214503000125}
#'
#' Hofner, B., Mayr, A., Robinzonov, N. and Schmid, M. (2014). Model-based
#' boosting in R: A hands-on tutorial using the R package mboost.
#' *Computational Statistics*, 29(1--2):3--35.
#' \doi{10.1007/s00180-012-0382-5}
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 4), n, 4)
#' y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
#' fit <- bgam(x, y, mstop = 50, autotune = FALSE)
#' plot(fit)
#' }
#' @export
plot.bgam <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2L, 2L))

  sf <- sort(x$selection_frequency, decreasing = TRUE)
  top_n <- min(4L, x$p)
  top_preds <- names(sf)[seq_len(top_n)]

  # Panel 1: Partial effect curve for top predictor (or first one)
  if (length(top_preds) >= 1L) {
    j_nm <- top_preds[1L]
    bl  <- x$base_learners[[j_nm]]
    bv  <- as.numeric(x$coefficients[[j_nm]])
    xgrid <- seq(bl$range[1L], bl$range[2L], length.out = 100L)
    if (length(bl$knots) == 2L && bl$K == 1L) {
      B_grid <- matrix(xgrid - mean(xgrid), ncol = 1L)
    } else {
      B_grid <- bgam_bspline_basis_cpp(xgrid, bl$knots, bl$degree)$B
    }
    fj_grid <- as.numeric(B_grid %*% bv)
    graphics::plot(xgrid, fj_grid, type = "l", col = "steelblue",
                   xlab = j_nm, ylab = paste0("f(", j_nm, ")"),
                   main = paste0("Partial effect: ", j_nm))
    graphics::abline(h = 0, lty = 2, col = "gray")
  }

  # Panel 2: Selection frequency bar chart
  bar_names <- names(sf)
  bar_vals  <- unname(sf)
  if (length(bar_names) > 20L) {
    bar_names <- bar_names[seq_len(20L)]
    bar_vals  <- bar_vals[seq_len(20L)]
  }
  graphics::barplot(bar_vals, names.arg = bar_names, las = 2L,
                    col = "steelblue", ylab = "Selection frequency",
                    main = "Predictor selection", cex.names = 0.7)

  # Panel 3: Training loss path
  iters  <- seq_along(x$loss_path)
  metric <- if (identical(x$family, "binomial")) "Deviance/n" else "MSE"
  graphics::plot(iters, x$loss_path, type = "l", col = "darkorange",
                 xlab = "Boosting iteration", ylab = metric,
                 main = "Training loss path")
  if (!is.null(x$cv))
    graphics::abline(v = x$cv$mstop_opt, lty = 2, col = "red")

  # Panel 4: Fitted vs observed / fitted by class
  if (identical(x$family, "gaussian")) {
    obs <- x$fitted.values + x$residuals
    graphics::plot(x$fitted.values, obs,
                   xlab = "Fitted", ylab = "Observed",
                   main = "Fitted vs Observed", pch = 20L, col = "gray40")
    graphics::abline(0, 1, lty = 2, col = "red")
  } else {
    # Recover response class from deviance residuals (sign gives class)
    ytr <- as.integer(x$residuals >= 0)
    graphics::boxplot(x$fitted.values ~ ytr,
                      xlab = "Observed class", ylab = "Fitted probability",
                      main = "Fitted prob. by class",
                      col = c("salmon", "steelblue"))
  }

  invisible(x)
}

# ----------------------------------------------------------------------------
#  plot.meep() helpers
# ----------------------------------------------------------------------------

# Out-of-fold R^2 over the rows where both response and prediction are finite.
# Returns NA when fewer than two usable rows or the response has zero variance.
.meep_r2 <- function(y, pred) {
  ok <- is.finite(y) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  yy <- y[ok]; pp <- pred[ok]
  ss_tot <- sum((yy - mean(yy))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - sum((yy - pp)^2) / ss_tot
}

# Tie-safe ROC curve + AUC for a 0/1 label vector scored by `score`.
# Returns list(fpr, tpr, auc); auc is NA (and the curve is degenerate) when
# either class is absent among the finite rows.
.meep_roc <- function(score, label) {
  ok <- is.finite(score) & is.finite(label)
  score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1)
  n0 <- sum(label == 0)
  if (n1 == 0L || n0 == 0L)
    return(list(fpr = c(0, 1), tpr = c(0, 1), auc = NA_real_))
  ord <- order(score, decreasing = TRUE)
  s <- label[ord]
  tpr <- cumsum(s) / sum(s)
  fpr <- cumsum(1 - s) / sum(1 - s)
  r <- rank(score)
  auc <- (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  list(fpr = c(0, fpr), tpr = c(0, tpr), auc = auc)
}

# Fixed palette recycled across learners (stack is always drawn in black).
.meep_palette <- function(n) {
  base <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7",
            "#E69F00", "#56B4E9", "#F0E442", "#999999")
  rep(base, length.out = max(1L, n))
}

# Names of learners "present" for a nuisance: columns with any finite entry.
.meep_present_learners <- function(mat) {
  if (is.null(mat) || !is.matrix(mat)) return(character(0))
  keep <- apply(mat, 2L, function(col) any(is.finite(col)))
  colnames(mat)[keep]
}

# Recover the response + stack prediction for one nuisance from a meep fit.
# Returns NULL when the nuisance is not present on the object.
.meep_nuisance_data <- function(x, nm) {
  if (identical(nm, "outcome")) {
    if (is.null(x$y_hat_oof)) return(NULL)
    resp <- x$y_resid + x$y_hat_oof
    return(list(resp = resp, stack = x$y_hat_oof,
                family = x$family, mat = x$oof_matrix$outcome))
  }
  if (identical(nm, "treatment")) {
    if (is.null(x$d_hat_oof) || is.null(x$treatment_family)) return(NULL)
    resp <- x$d_resid + x$d_hat_oof
    return(list(resp = resp, stack = x$d_hat_oof,
                family = x$treatment_family, mat = x$oof_matrix$treatment))
  }
  if (identical(nm, "mu0") || identical(nm, "mu1")) {
    stack <- if (identical(nm, "mu0")) x$mu0_hat_oof else x$mu1_hat_oof
    if (is.null(stack)) return(NULL)
    # mu0 / mu1 responses are the outcome response on the arm rows.
    resp <- x$y_resid + x$y_hat_oof
    return(list(resp = resp, stack = stack,
                family = x$family, mat = x$oof_matrix[[nm]]))
  }
  NULL
}

# Draw an empty panel carrying an explanatory message.
.meep_empty_panel <- function(main, msg) {
  graphics::plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
                 main = main)
  graphics::text(0, 0, msg, col = "gray40")
}

# Draw the ROC panel for a binomial nuisance.
.meep_roc_panel <- function(nd, nm) {
  resp <- nd$resp; stack <- nd$stack
  present <- .meep_present_learners(nd$mat)
  stack_roc <- .meep_roc(stack, resp)
  curves <- lapply(present, function(l) .meep_roc(nd$mat[, l], resp))
  drawable <- vapply(curves, function(r) !is.na(r$auc), logical(1))
  if (is.na(stack_roc$auc) && !any(drawable)) {
    .meep_empty_panel(paste0("ROC: ", nm),
                      "ROC undefined (single-class response)")
    return(invisible())
  }
  graphics::plot(c(0, 1), c(0, 1), type = "n",
                 xlab = "False positive rate", ylab = "True positive rate",
                 main = paste0("ROC: ", nm))
  graphics::abline(0, 1, lty = 3, col = "gray60")
  cols <- .meep_palette(length(present))
  leg_txt <- character(0); leg_col <- character(0); leg_lwd <- numeric(0)
  for (i in seq_along(present)) {
    r <- curves[[i]]
    if (is.na(r$auc)) next
    graphics::lines(r$fpr, r$tpr, col = cols[i], lwd = 1)
    leg_txt <- c(leg_txt, sprintf("%s (%.2f)", present[i], r$auc))
    leg_col <- c(leg_col, cols[i]); leg_lwd <- c(leg_lwd, 1)
  }
  if (!is.na(stack_roc$auc)) {
    graphics::lines(stack_roc$fpr, stack_roc$tpr, col = "black", lwd = 2.5)
    leg_txt <- c(leg_txt, sprintf("stack (%.2f)", stack_roc$auc))
    leg_col <- c(leg_col, "black"); leg_lwd <- c(leg_lwd, 2.5)
  }
  if (length(leg_txt))
    graphics::legend("bottomright", legend = leg_txt, col = leg_col,
                     lwd = leg_lwd, bty = "n", cex = 0.8)
  invisible()
}

# Draw the two gaussian panels (R^2 barplot + observed-vs-fitted scatter).
.meep_gaussian_panels <- function(nd, nm) {
  resp <- nd$resp; stack <- nd$stack
  present <- .meep_present_learners(nd$mat)
  # Panel 1: OOF R^2 barplot (per learner + stack).
  r2_learner <- vapply(present, function(l) .meep_r2(resp, nd$mat[, l]),
                       numeric(1))
  r2_stack <- .meep_r2(resp, stack)
  heights <- c(r2_learner, stack = r2_stack)
  heights[!is.finite(heights)] <- 0
  cols <- c(.meep_palette(length(present)), "black")
  graphics::barplot(heights, names.arg = names(heights), las = 2L,
                    col = cols, ylab = "OOF R^2", cex.names = 0.8,
                    main = paste0("OOF R^2: ", nm))
  graphics::abline(h = 0, col = "gray70", lwd = 1)
  # Panel 2: observed vs OOF stack prediction.
  ok <- is.finite(resp) & is.finite(stack)
  if (!any(ok)) {
    .meep_empty_panel(paste0("Observed vs OOF fit: ", nm),
                      "no finite rows")
    return(invisible())
  }
  graphics::plot(stack[ok], resp[ok], pch = 20L, col = "gray40",
                 xlab = "OOF prediction", ylab = "Response",
                 main = paste0("Observed vs OOF fit: ", nm))
  graphics::abline(0, 1, lty = 2, col = "red")
  if (is.finite(r2_stack))
    graphics::legend("topleft", legend = sprintf("stack R^2 = %.2f", r2_stack),
                     bty = "n", cex = 0.85)
  invisible()
}

#' Diagnostic plots for a `meep` cross-fitted ensemble
#'
#' Visualises out-of-fold (OOF) ensemble performance per nuisance.  For each
#' selected nuisance the panel layout is keyed to that nuisance's family:
#' \itemize{
#'   \item \strong{binomial} -- one ROC panel overlaying each present learner's
#'     OOF ROC curve (thin, coloured) plus the combined "stack" curve (thick,
#'     black), with a 45-degree reference line and a legend reporting each
#'     curve's AUC.
#'   \item \strong{gaussian} -- two panels: an OOF \eqn{R^2} bar chart (one bar
#'     per present learner plus a "stack" bar) and an observed-vs-OOF-prediction
#'     scatter for the stack with a 45-degree line and the stack \eqn{R^2}
#'     annotated.
#' }
#'
#' Panels are laid out on a near-square grid sized to the total panel count
#' (binomial contributes 1 panel, gaussian 2).  The function works in headless
#' (non-interactive) mode without error; degenerate nuisances (single-class
#' response, zero-variance outcome) are skipped gracefully rather than erroring.
#'
#' @param x An object of class `"meep"`, as returned by [meep()].
#' @param which Character vector of nuisances to draw; any subset of
#'   `c("outcome", "treatment", "mu0", "mu1")`.  Nuisances not present on the
#'   object are silently dropped.  Defaults to `c("outcome", "treatment")`.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @seealso [meep()]
#' @export
plot.meep <- function(x, which = c("outcome", "treatment"), ...) {
  if (!inherits(x, "meep"))
    stop("plot.meep: 'x' must be a 'meep' object.")
  allowed <- c("outcome", "treatment", "mu0", "mu1")
  which <- intersect(which, allowed)
  if (length(which) == 0L)
    stop("plot.meep: `which` must name at least one of ",
         paste(allowed, collapse = ", "), ".")

  # Resolve which selected nuisances are actually present.
  nds <- list()
  for (nm in which) {
    nd <- .meep_nuisance_data(x, nm)
    if (!is.null(nd)) nds[[nm]] <- nd
  }
  if (length(nds) == 0L)
    stop("plot.meep: none of the requested nuisances (",
         paste(which, collapse = ", "), ") are present on this fit.")

  # Total panel count: binomial = 1, gaussian = 2.
  n_panels <- sum(vapply(nds, function(nd)
    if (identical(nd$family, "binomial")) 1L else 2L, integer(1)))
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))

  for (nm in names(nds)) {
    nd <- nds[[nm]]
    if (identical(nd$family, "binomial")) {
      .meep_roc_panel(nd, nm)
    } else {
      .meep_gaussian_panels(nd, nm)
    }
  }

  invisible(x)
}
