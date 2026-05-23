# BUG-015 regression test (v0.0.0.9055).
#
# predict.ares() silently dropped newdata rows with NA in numeric columns
# when the fit had a $factor_info expansion. The factor_info branch
# called model.matrix(~ ., data = newdata) with the global na.action
# (default na.omit), which removed offending rows -- the returned vector
# length was less than nrow(newdata), violating the predict()-length
# contract.
#
# Fix: build the model frame with na.action = na.pass, then model.matrix,
# so NAs in numeric columns flow through to the design and the
# downstream median-impute / na-row-mask path handles them.

test_that("predict.ares with NA in numeric newdata returns nrow(newdata)", {
  set.seed(1)
  n <- 50
  df <- data.frame(
    x = rnorm(n),
    g = factor(sample(c("a", "b"), n, replace = TRUE)),
    y = rnorm(n)
  )
  fit <- ares(x = df[, c("x", "g")], y = df$y)
  new_df <- data.frame(x = c(NA, 0.5, NA),
                       g = factor(c("a", "b", "a"), levels = c("a", "b")))
  # warning is expected for the NA imputation path; length must be 3.
  out <- suppressWarnings(predict(fit, new_df))
  expect_equal(length(out), 3L)
})

test_that("predict.ares formula path with NA in numeric newdata is nrow-preserving", {
  set.seed(2)
  n <- 50
  df <- data.frame(
    x = rnorm(n),
    g = factor(sample(c("a", "b"), n, replace = TRUE)),
    y = rnorm(n)
  )
  fit <- ares(y ~ x + g, data = df)
  new_df <- data.frame(x = c(NA, 0.5, NA),
                       g = factor(c("a", "b", "a"), levels = c("a", "b")))
  out <- suppressWarnings(predict(fit, new_df))
  expect_equal(length(out), 3L)
})
