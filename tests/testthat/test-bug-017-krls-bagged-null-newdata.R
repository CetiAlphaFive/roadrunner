# BUG-017 regression test (v0.0.0.9055).
#
# predict(krls_bagged_fit, newdata = NULL) returned $fitted -- the
# central (single) fit -- instead of the bag-mean across the bootstrap
# replicates. The non-NULL branch already does the bagging, so the
# NULL-newdata fast-path was inconsistent with the predict-on-newdata
# behaviour. CLAUDE.md and predict.ares both promise bag-mean from
# NULL-newdata.
#
# Fix: in the NULL-newdata branch, when $boot$replicates exists, route
# through the non-NULL path with newdata = object$X.

test_that("predict(krls_bagged) with NULL newdata equals predict on training X", {
  set.seed(7)
  n <- 30
  X <- matrix(rnorm(n * 2), n, 2)
  y <- sin(X[, 1]) + 0.2 * rnorm(n)
  fit <- krls(X, y, n.boot = 4, derivative = FALSE, vcov = FALSE)
  p_null <- predict(fit)
  p_full <- predict(fit, X)
  # Both objects must be the same shape (list with $fit) and numerically
  # equal -- the central fit alone would NOT match the bag mean.
  expect_equal(as.numeric(p_null$fit), as.numeric(p_full$fit),
               tolerance = 1e-12)
})
