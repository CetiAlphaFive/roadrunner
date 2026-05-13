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
