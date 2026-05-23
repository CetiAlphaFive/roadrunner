# BUG-016 regression test (v0.0.0.9055).
#
# Weighted krls() returned posterior variance ~= 0 from
# predict(..., type = "variance"). The fit stored (V, dvals) from the
# WEIGHTED K_w = D K D (D = diag(sqrt(w))), but the predict-side helper
# built K_star from the UNWEIGHTED kernel. The two bases were
# inconsistent: for the Gaussian kernel where diag(K**) = 1 the
# quad-form blew past 1 and the floor-at-zero hid the bug as v == 0
# everywhere.
#
# Fix: in .krls_predict_variance, when training weights are present,
# scale K_star columnwise by sqrt(w_tr) before contracting against
# (V, dvals). The math identity is
#   V[f*] = K**(unwt) - K* D (K_w + lam I)^-1 D K*'.

test_that("predict(krls weighted fit, type='variance') is non-zero", {
  set.seed(11)
  n <- 30
  X <- matrix(rnorm(n * 2), n, 2)
  y <- sin(X[, 1]) + 0.2 * rnorm(n)
  w <- runif(n) + 0.1
  fit <- krls(X, y, weights = w, derivative = FALSE, vcov = TRUE)
  v <- predict(fit, matrix(rnorm(6), 3, 2), type = "variance")
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0))
  # The bug rendered v exactly == 0 across the whole vector under the
  # floor-at-zero clamp. After the fix every entry should be strictly
  # positive for any reasonable newdata row.
  expect_true(all(v > 0))
})

test_that("unweighted krls posterior variance is unchanged by the fix", {
  set.seed(12)
  n <- 25
  X <- matrix(rnorm(n * 2), n, 2)
  y <- sin(X[, 1]) + 0.2 * rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = TRUE)
  v <- predict(fit, matrix(rnorm(6), 3, 2), type = "variance")
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0))
})
