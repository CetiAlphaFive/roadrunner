test_that("plda fit is byte-identical across thread counts", {
  set.seed(8)
  x <- matrix(rnorm(120 * 15), 120, 15)
  y <- factor(rep(1:3, length.out = 120))
  f1 <- plda(x, y, K = 2, lambda = 0.07, autotune = FALSE, nthreads = 1L)
  f4 <- plda(x, y, K = 2, lambda = 0.07, autotune = FALSE, nthreads = 4L)
  expect_identical(f1$discrim, f4$discrim)
})

test_that("plda autotune is byte-identical across thread counts", {
  set.seed(9)
  x <- matrix(rnorm(150 * 12), 150, 12)
  y <- factor(rep(1:3, length.out = 150))
  f1 <- plda(x, y, nthreads = 1L)
  f4 <- plda(x, y, nthreads = 4L)
  expect_identical(f1$discrim, f4$discrim)
  expect_identical(f1$lambda, f4$lambda)
})
