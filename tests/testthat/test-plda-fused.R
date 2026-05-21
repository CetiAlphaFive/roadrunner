test_that("plda_tv1d_cpp solves the 1-D fused-lasso signal approximator", {
  u <- c(1, 1, 1, 5, 5, 5)
  got <- roadrunner:::plda_tv1d_cpp(u, 0.0)
  expect_equal(got, u, tolerance = 1e-10)
  big <- roadrunner:::plda_tv1d_cpp(u, 100)
  expect_equal(big, rep(mean(u), 6), tolerance = 1e-6)
})

test_that("plda fused penalty matches penalizedLDA type='ordered'", {
  skip_if_not_installed("penalizedLDA")
  set.seed(99)
  n <- 80; p <- 30
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:2, length.out = n)
  x[y == 1, 10:20] <- x[y == 1, 10:20] + 1.2
  ref <- penalizedLDA::PenalizedLDA(x, y, lambda = 0.05, K = 1,
                                    type = "ordered", lambda2 = 0.1)
  fit <- plda(x, factor(y), K = 1, lambda = 0.05, penalty = "fused",
              lambda2 = 0.1, autotune = FALSE)
  a <- fit$discrim[, 1]; b <- ref$discrim[, 1]
  if (sum((a - b)^2) > sum((a + b)^2)) a <- -a
  expect_equal(a, b, tolerance = 1e-3)
})
