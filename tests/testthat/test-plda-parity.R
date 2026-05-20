test_that("plda_fit_cpp K=1 L1 matches penalizedLDA::PenalizedLDA", {
  skip_if_not_installed("penalizedLDA")
  set.seed(42)
  n <- 60; p <- 20
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:3, length.out = n)
  x[y == 1, 1:3] <- x[y == 1, 1:3] + 2
  ref <- penalizedLDA::PenalizedLDA(x, y, lambda = 0.1, K = 1, type = "standard")
  got <- roadrunner:::plda_fit_cpp(x, as.integer(y), 3L, 1L, 0.1, 0.0,
                                   0L, 100L, 1e-6)
  d_got <- got$discrim[, 1]; d_ref <- ref$discrim[, 1]
  if (sum((d_got - d_ref)^2) > sum((d_got + d_ref)^2)) d_got <- -d_got
  expect_equal(d_got, d_ref, tolerance = 1e-4)
})
