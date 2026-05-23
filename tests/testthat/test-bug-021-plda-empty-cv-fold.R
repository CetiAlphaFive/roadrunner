# BUG-021 regression test (v0.0.0.9055).
#
# plda() with nfold > min(class_count) silently produced biased CV
# error. The stratified fold builder placed each class into folds
# 1..length(class) only, so when min(class_count) < nfold some folds got
# zero test observations; those folds short-circuited with err_slot = 0
# in the C++ harness, biasing the average CV error toward 0 and
# selecting the wrong lambda.
#
# Fix: reject upfront when nfold > min(table(y)) with a clear message
# recommending the largest valid nfold.

test_that("plda() errors when nfold exceeds the smallest class count", {
  set.seed(1)
  x <- matrix(rnorm(4 * 5), 4, 5)
  y <- factor(c(1, 1, 2, 2))   # 2 obs/class
  expect_error(plda(x, y, nfold = 4L), "nfold")
})

test_that("plda() still fits when nfold == min(class_count)", {
  set.seed(2)
  # 8 obs/class with nfold = 4 -> each fold gets 2 obs/class.
  n_per <- 8L
  x <- rbind(matrix(rnorm(n_per * 5), n_per, 5) + 0.5,
             matrix(rnorm(n_per * 5), n_per, 5) - 0.5)
  y <- factor(rep(c(1, 2), each = n_per))
  expect_no_error(plda(x, y, nfold = 4L))
})
