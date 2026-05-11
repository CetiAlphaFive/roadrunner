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
              glm = object$glm)
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
  cat("\n")
  print(x$terms, row.names = FALSE)
  invisible(x)
}

#' Plot method for `ares` fits
#'
#' Diagnostic plot for a fitted `ares` model. Currently produces a
#' residuals-vs-fitted scatter; additional diagnostics may be added in
#' the future.
#'
#' @param x An object of class `"ares"`.
#' @param which Which diagnostic to draw. Only `1` (residuals vs
#'   fitted) is implemented.
#' @param ... Further graphical parameters passed to `plot()`.
#' @return Invisibly returns `x`.
#' @examples
#' \dontrun{
#'   fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#'   plot(fit)
#' }
#' @export
plot.ares <- function(x, which = 1L, ...) {
  graphics::plot(x$fitted.values, x$residuals,
                 xlab = "Fitted", ylab = "Residual",
                 main = "ares: residuals vs fitted", ...)
  graphics::abline(h = 0, lty = 2)
  invisible(x)
}
