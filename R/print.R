# S3 methods: print, summary, print.summary, plot.
#
# All four are thin formatters over fields populated by ares.default():
#   $call, $dirs, $selected.terms, $rss, $gcv, $degree, $penalty,
#   $nthreads, $coefficients, $fitted.values, $residuals.

#' Print method for `ares` fits
#' @param x an `ares` object
#' @param digits significant digits for numeric output
#' @param ... unused
#' @return Invisibly returns `x`.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' print(fit)
#' @export
print.ares <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMARS fit:\n")
  cat("  Selected terms:", length(x$selected.terms), "of", nrow(x$dirs),
      "forward-pass terms\n")
  cat("  RSS:", format(x$rss, digits = digits), "\n")
  cat("  GCV:", format(x$gcv, digits = digits), "\n")
  cat("  Degree:", x$degree, "  Penalty:", x$penalty, "  nthreads:", x$nthreads, "\n")
  invisible(x)
}

#' Summary method for `ares` fits
#' @param object an `ares` object
#' @param ... unused
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
  out <- list(call = object$call, terms = tab,
              rss = object$rss, gcv = object$gcv,
              n_terms = length(object$selected.terms),
              n_forward = nrow(object$dirs),
              degree = object$degree, penalty = object$penalty)
  class(out) <- "summary.ares"
  out
}

#' @rdname summary.ares
#' @param x a `summary.ares` object
#' @export
print.summary.ares <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nMARS summary:\n")
  cat("  Selected:", x$n_terms, "of", x$n_forward, "forward-pass terms\n")
  cat("  RSS =", format(x$rss), "  GCV =", format(x$gcv),
      "  degree =", x$degree, "\n\n")
  print(x$terms, row.names = FALSE)
  invisible(x)
}

#' Plot method for `ares` fits (residuals vs fitted)
#' @param x an `ares` object
#' @param which integer in 1:1 — currently only residuals-vs-fitted is supported
#' @param ... passed to `plot()`
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
