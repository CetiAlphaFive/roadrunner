# BUG-014 regression test (v0.0.0.9055).
#
# ares(pmethod="cv") crashed with "subscript out of bounds" when a small
# fold's forward pass collapsed to an intercept-only model (M_full == 1).
# The C++ engine only populates $path.subsets / $path.coefs when
# pmethod == 0 && M > 1, so those lists are empty for an M=1 fold.
# The R-side CV loop indexed `fit$path.subsets[[s]]` without checking
# list length, throwing "subscript out of bounds".
#
# Fix: guard the index against length(fit$path.subsets) /
# length(fit$path.coefs) before extracting.

test_that("ares(pmethod='cv') survives intercept-only fold collapse", {
  set.seed(1)
  n <- 10
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  expect_no_error(
    fit <- ares(x, y, pmethod = "cv", nfold = 5, seed.cv = 1, nthreads = 1)
  )
  expect_s3_class(fit, "ares")
  # And predict() returns a length-n vector.
  expect_equal(length(predict(fit, x)), n)
})
