test_that("plda_tv1d_cpp solves the 1-D fused-lasso signal approximator", {
  u <- c(1, 1, 1, 5, 5, 5)
  got <- roadrunner:::plda_tv1d_cpp(u, 0.0)
  expect_equal(got, u, tolerance = 1e-10)
  big <- roadrunner:::plda_tv1d_cpp(u, 100)
  expect_equal(big, rep(mean(u), 6), tolerance = 1e-6)
})
