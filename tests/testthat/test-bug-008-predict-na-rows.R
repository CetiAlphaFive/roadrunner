# BUG-008 regression test (v0.0.0.9032).
#
# predict() used to return finite WRONG values for newdata rows with any
# NA when training used na.action = "omit". The warning at predict.R:184
# promised "the affected rows will return NA predictions" but the code
# never set those rows to NA. NaN passed to the C++ hinge `(x - t) > 0`
# evaluates to FALSE (NaN comparisons are always false), so each hinge
# involving an NA variable collapsed to 0 and the row got a deterministic
# but wrong prediction (intercept plus the non-NA contributions).
#
# Fix: detect rows with any NA in xnew BEFORE mars_basis_cpp, zero-fill
# the NaN cells so the C++ basis pass sees finite input, then re-impose
# NA on those rows post linear-predictor compute. Apply on yhat, on bag
# SE, and on the interval = "pint" matrix path.

test_that("predict() returns NA for newdata rows with NA (matrix interface, omit)", {
  set.seed(2)
  x <- matrix(rnorm(180), 60, 3)
  y <- x[, 1] + rnorm(60) * 0.5
  fit <- ares(x, y, na.action = "omit", nthreads = 2)
  x_pred <- x
  x_pred[1:5, 1] <- NA_real_
  expect_warning(p <- predict(fit, x_pred),
                 "Missing values in newdata")
  expect_true(all(is.na(p[1:5])))
  expect_true(all(!is.na(p[6:60])))
  expect_equal(length(p), nrow(x_pred))
})

test_that("predict() returns NA for newdata rows with NA spread across columns", {
  set.seed(3)
  x <- matrix(rnorm(240), 80, 3)
  y <- x[, 1] - x[, 2] + rnorm(80) * 0.3
  fit <- ares(x, y, na.action = "omit", nthreads = 2)
  x_pred <- x
  x_pred[10, 1] <- NA
  x_pred[20, 2] <- NA
  x_pred[30, 3] <- NA
  suppressWarnings(p <- predict(fit, x_pred))
  expect_true(all(is.na(p[c(10, 20, 30)])))
  expect_true(all(!is.na(p[-c(10, 20, 30)])))
})

test_that("predict() returns NA on PI matrix for NA rows (omit + varmod)", {
  set.seed(4)
  x <- matrix(rnorm(180), 60, 3)
  y <- x[, 1] + rnorm(60) * 0.5
  fit <- ares(x, y, na.action = "omit", varmod = "const", nthreads = 2)
  x_pred <- x
  x_pred[1:3, 1] <- NA
  suppressWarnings(pm <- predict(fit, x_pred, interval = "pint"))
  expect_true(is.matrix(pm))
  expect_equal(ncol(pm), 3L)
  expect_true(all(is.na(pm[1:3, ])))
  expect_true(all(!is.na(pm[4:60, ])))
})

test_that("predict() na.action='impute' path still imputes (no regression)", {
  set.seed(5)
  x <- matrix(rnorm(120), 40, 3)
  # Introduce a training NA so the fit stashes na.medians (which is what
  # triggers the impute branch in predict()).
  x[1, 1] <- NA_real_
  y <- rowSums(x, na.rm = TRUE) * 0.3 + rnorm(40) * 0.5
  fit <- ares(x, y, na.action = "impute", nthreads = 2)
  expect_false(is.null(fit$na.medians))
  x_pred <- x
  x_pred[2:4, 1] <- NA
  expect_warning(p <- predict(fit, x_pred), "median-imputed")
  expect_true(all(!is.na(p)))
})
