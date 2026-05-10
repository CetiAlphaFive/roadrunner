# Pure additive sum-of-hinges DGP

#' @noRd
dgp_additive <- function(n, p = 10, snr = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  a <- c(2, 1.5, 1, 0.5, -1)
  active_p <- min(5, p)
  f <- numeric(n)
  for (j in seq_len(active_p)) f <- f + a[j] * pmax(0, x[, j] - 0.5)
  sigma <- sqrt(stats::var(f) / snr)
  y <- f + stats::rnorm(n, sd = sigma)
  list(x = x, y = y, f = f, sigma = sigma)
}
