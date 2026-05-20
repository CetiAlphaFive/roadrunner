test_that("plda_wcsd_cpp returns within-class sd with divisor n", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  g <- rep(1:2, each = 10)
  got <- roadrunner:::plda_wcsd_cpp(x, g, 2L)
  expect_length(got, 2L)
  ref <- vapply(1:2, function(j) {
    ss <- sum(tapply(seq_len(20), g, function(ix) sum((x[ix, j] - mean(x[ix, j]))^2)))
    sqrt(ss / 20)
  }, numeric(1))
  expect_equal(got, ref, tolerance = 1e-10)
})
