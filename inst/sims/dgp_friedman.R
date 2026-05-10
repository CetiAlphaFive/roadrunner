# Friedman 1991 / 1992 canonical DGPs for MARS evaluation.
# Reproducibility: each generator takes a `seed` argument.

#' Friedman-1 DGP
#'
#' y = 10 sin(pi x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5 + eps
#' x ~ U(0, 1)^p, with p >= 5 (extra columns are nuisance).
#'
#' @param n integer; number of observations
#' @param p integer; total number of predictors (>= 5)
#' @param snr numeric; signal-to-noise ratio (Var(f)/sigma^2)
#' @param seed integer or NULL
#' @return list(x, y, f, sigma)
#' @noRd
dgp_friedman1 <- function(n, p = 10, snr = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  f <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5]
  sigma <- sqrt(stats::var(f) / snr)
  y <- f + stats::rnorm(n, sd = sigma)
  list(x = x, y = y, f = f, sigma = sigma)
}

#' Friedman-2 DGP (atan-style)
#' @noRd
dgp_friedman2 <- function(n, p = 5, snr = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::runif(n, 0, 100)
  x2 <- stats::runif(n, 40 * pi, 560 * pi)
  x3 <- stats::runif(n, 0, 1)
  x4 <- stats::runif(n, 1, 11)
  f <- sqrt(x1^2 + (x2 * x3 - 1 / (x2 * x4))^2)
  X <- cbind(x1, x2, x3, x4)
  if (p > 4) X <- cbind(X, matrix(stats::runif(n * (p - 4)), n, p - 4))
  sigma <- sqrt(stats::var(f) / snr)
  y <- f + stats::rnorm(n, sd = sigma)
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  list(x = X, y = y, f = f, sigma = sigma)
}

#' Friedman-3 DGP
#' @noRd
dgp_friedman3 <- function(n, p = 5, snr = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::runif(n, 0, 100)
  x2 <- stats::runif(n, 40 * pi, 560 * pi)
  x3 <- stats::runif(n, 0, 1)
  x4 <- stats::runif(n, 1, 11)
  f <- atan2(x2 * x3 - 1 / (x2 * x4), x1)
  X <- cbind(x1, x2, x3, x4)
  if (p > 4) X <- cbind(X, matrix(stats::runif(n * (p - 4)), n, p - 4))
  sigma <- sqrt(stats::var(f) / snr)
  y <- f + stats::rnorm(n, sd = sigma)
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  list(x = X, y = y, f = f, sigma = sigma)
}
