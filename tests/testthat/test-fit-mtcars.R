# Mtcars fit-quality test. ares is no longer required to match earth;
# we only check that the fit is not pathological on a small classic data set.

test_that("ares fits mtcars with sensible RSS / R^2", {
  skip_on_cran()
  x <- as.matrix(mtcars[, -1]); y <- mtcars$mpg
  fa <- ares(x, y, nthreads = 2)
  pa <- predict(fa, x)
  ss_tot <- sum((y - mean(y))^2)
  ss_res <- sum((y - pa)^2)
  r2 <- 1 - ss_res / ss_tot
  # Mtcars n=32 is tiny and rank-deficient; reasonable training-set R^2 ~ 0.85.
  expect_gt(r2, 0.85)
  # GCV finite and positive.
  expect_true(is.finite(fa$gcv) && fa$gcv > 0)
})
