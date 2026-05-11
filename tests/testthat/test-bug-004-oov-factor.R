# BUG-004 regression test (v0.0.0.9029).
#
# OOV factor levels in newdata used to become NA, then model.matrix()
# silently dropped those rows, so predict() returned fewer rows than
# nrow(newdata). predict.ares() now detects OOV up front and stop()s with
# a clear, actionable error message.

test_that("predict errors loudly on OOV factor levels in newdata (factor col)", {
  set.seed(1)
  df <- data.frame(num = stats::runif(50),
                   cat = factor(sample(c("aaa", "bb-bb"), 50, replace = TRUE)))
  yy <- df$num + (df$cat == "bb-bb") + stats::rnorm(50, sd = 0.3)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$cat <- as.character(df2$cat)
  df2$cat[1:5] <- "OOV"
  expect_error(predict(fit, df2),
               "out-of-vocabulary factor level")
})

test_that("predict errors loudly on OOV factor levels in newdata (character col)", {
  set.seed(2)
  df <- data.frame(num = stats::runif(40),
                   cat = sample(c("A", "B", "C"), 40, replace = TRUE),
                   stringsAsFactors = FALSE)
  yy <- df$num + as.integer(df$cat == "B") + stats::rnorm(40, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$cat[c(2, 7, 11)] <- "ZZZ"
  expect_error(predict(fit, df2),
               "out-of-vocabulary factor level")
})

test_that("in-vocabulary newdata still returns nrow(newdata) predictions", {
  set.seed(3)
  df <- data.frame(num = stats::runif(40),
                   cat = factor(sample(c("a", "b"), 40, replace = TRUE)))
  yy <- df$num + (df$cat == "b") + stats::rnorm(40, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  p <- predict(fit, df)
  expect_equal(length(p), nrow(df))
})

test_that("OOV error message names the column and offending level", {
  set.seed(4)
  df <- data.frame(num = stats::runif(30),
                   cat = factor(sample(c("p", "q"), 30, replace = TRUE)))
  yy <- df$num + (df$cat == "q") + stats::rnorm(30, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$cat <- as.character(df2$cat)
  df2$cat[1] <- "BADLEVEL"
  msg <- tryCatch(predict(fit, df2),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "column cat")
  expect_match(msg, "BADLEVEL")
})
