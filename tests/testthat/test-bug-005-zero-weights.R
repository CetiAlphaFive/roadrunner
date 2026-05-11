# BUG-005 regression test (v0.0.0.9029).
#
# Zero weights used to be permitted, but the C++ engine's GCV denominator
# and varmod df assume sum(w) = n_effective via the mean(w) = 1
# normalisation. With zero-weight rows present, sum(w) = n while
# n_effective = sum(w > 0) < n, so the GCV penalty was implicitly
# under-estimated and the model over-fit. We now reject any zero (or
# negative) weight up front and tell the user to drop those rows
# instead.

test_that("zero weights are rejected with a clear error", {
  set.seed(120)
  n <- 100; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  y <- 2 * x[, 1] + stats::rnorm(n, sd = 0.3)
  w <- runif(n, 0.5, 2)
  w_one_zero <- w
  w_one_zero[5] <- 0
  expect_error(ares(x, y, weights = w_one_zero, nthreads = 2),
               "strictly positive")
})

test_that("all-zero weights are rejected with a clear error", {
  set.seed(121)
  n <- 80; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  y <- x[, 1] + stats::rnorm(n, sd = 0.3)
  expect_error(ares(x, y, weights = rep(0, n), nthreads = 2),
               "strictly positive")
})

test_that("dropping zero-weight rows yields a valid fit (lossless workaround)", {
  set.seed(122)
  n <- 200; p <- 4
  x <- matrix(stats::runif(n * p), n, p)
  y <- 3 * x[, 1] + sin(2 * pi * x[, 2]) + stats::rnorm(n, sd = 0.4)
  w <- runif(n, 0.5, 2)
  w[1:50] <- 0
  expect_error(ares(x, y, weights = w, nthreads = 2),
               "strictly positive")
  # User-side workaround: drop the rows.
  keep <- w > 0
  fdr <- ares(x[keep, , drop = FALSE], y[keep],
              weights = w[keep], nthreads = 2)
  expect_true(is.finite(fdr$rss))
  expect_true(length(fdr$selected.terms) >= 1)
})

test_that("error message points the user to the correct workaround", {
  set.seed(123)
  n <- 60; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  y <- x[, 1] + stats::rnorm(n, sd = 0.3)
  w <- rep(1, n); w[1] <- 0
  msg <- tryCatch(ares(x, y, weights = w, nthreads = 2),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "drop", fixed = TRUE)
  expect_match(msg, "GCV")
})
