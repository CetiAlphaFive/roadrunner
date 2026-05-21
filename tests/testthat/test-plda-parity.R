test_that("plda_fit_cpp K=1 L1 matches penalizedLDA::PenalizedLDA", {
  skip_if_not_installed("penalizedLDA")
  PenalizedLDA <- getExportedValue("penalizedLDA", "PenalizedLDA")
  set.seed(42)
  n <- 60; p <- 20
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:3, length.out = n)
  x[y == 1, 1:3] <- x[y == 1, 1:3] + 2
  ref <- PenalizedLDA(x, y, lambda = 0.1, K = 1, type = "standard")
  got <- roadrunner:::plda_fit_cpp(x, as.integer(y), 3L, 1L, 0.1, 0.0,
                                   0L, 100L, 1e-6)
  d_got <- got$discrim[, 1]; d_ref <- ref$discrim[, 1]
  if (sum((d_got - d_ref)^2) > sum((d_got + d_ref)^2)) d_got <- -d_got
  expect_equal(d_got, d_ref, tolerance = 1e-4)
})

test_that("plda_fit_cpp K=2 L1 matches penalizedLDA discrim matrix", {
  skip_if_not_installed("penalizedLDA")
  PenalizedLDA <- getExportedValue("penalizedLDA", "PenalizedLDA")
  set.seed(7)
  n <- 90; p <- 25
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:3, length.out = n)
  x[y == 1, 1:4] <- x[y == 1, 1:4] + 1.5
  x[y == 2, 5:8] <- x[y == 2, 5:8] - 1.5
  ref <- PenalizedLDA(x, y, lambda = 0.08, K = 2, type = "standard")
  got <- roadrunner:::plda_fit_cpp(x, as.integer(y), 3L, 2L, 0.08, 0.0, 0L, 200L, 1e-7)
  for (k in 1:2) {
    a <- got$discrim[, k]; b <- ref$discrim[, k]
    if (sum((a - b)^2) > sum((a + b)^2)) a <- -a
    expect_equal(a, b, tolerance = 1e-3, info = paste("discriminant", k))
  }
})
