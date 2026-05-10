# Two-way interaction (degree-2) DGP

#' @noRd
dgp_interaction <- function(n, p = 10, snr = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  f <- 5 * pmax(0, x[, 1] - 0.5) * pmax(0, x[, 2] - 0.3) +
       3 * pmax(0, x[, 3] - 0.6)
  sigma <- sqrt(stats::var(f) / snr)
  y <- f + stats::rnorm(n, sd = sigma)
  list(x = x, y = y, f = f, sigma = sigma)
}
