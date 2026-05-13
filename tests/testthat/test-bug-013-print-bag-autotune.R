# BUG-013 regression test (v0.0.0.9032).
#
# print.ares() and summary.ares() used to be silent about bagging
# (n.boot > 0) and autotune state. A bagged or autotuned fit printed
# identically to a plain fit; users had to inspect the object structure
# manually to see the bag size or the autotune winner.
#
# Fix: add a one-line block to each printer when $boot or $autotune is
# present. The Call: line is left as-is.

skip_if_quick <- function() {
  if (!nzchar(Sys.getenv("ARES_FULL_TESTS"))) skip("ARES_FULL_TESTS not set")
}

test_that("print.ares mentions Bagging when n.boot > 0", {
  set.seed(1)
  x <- matrix(rnorm(60 * 3), 60, 3)
  y <- x[, 1] + rnorm(60, sd = 0.3)
  fit <- ares(x, y, n.boot = 3, nthreads = 2, seed.cv = 1)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Bagging", out)))
  expect_true(any(grepl("n\\.boot = 3", out)))
})

test_that("print.ares mentions Autotune when autotune = TRUE", {
  skip_if_quick()
  set.seed(2)
  x <- matrix(rnorm(120 * 3), 120, 3)
  y <- x[, 1] + rnorm(120, sd = 0.3)
  fit <- ares(x, y, autotune = TRUE, nthreads = 2, seed.cv = 2)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Autotune", out)))
  expect_true(any(grepl("degree=", out)))
})

test_that("summary.ares mentions Bagging when n.boot > 0", {
  set.seed(3)
  x <- matrix(rnorm(60 * 3), 60, 3)
  y <- x[, 1] + rnorm(60, sd = 0.3)
  fit <- ares(x, y, n.boot = 3, nthreads = 2, seed.cv = 3)
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Bagging", out)))
})

test_that("print.ares is silent about Bagging/Autotune for plain fit", {
  set.seed(4)
  x <- matrix(rnorm(60 * 3), 60, 3)
  y <- x[, 1] + rnorm(60, sd = 0.3)
  fit <- ares(x, y, nthreads = 2)
  out <- capture.output(print(fit))
  expect_false(any(grepl("Bagging", out)))
  expect_false(any(grepl("Autotune", out)))
})
