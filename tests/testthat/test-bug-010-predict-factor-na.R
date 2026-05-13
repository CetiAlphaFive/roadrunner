# BUG-010 regression test (v0.0.0.9032).
#
# Sister to BUG-004 (OOV factor levels). NA values in factor or character
# columns of newdata used to fall through the OOV detector (which only
# checked !is.na(.)), then model.matrix(~ ., newdata)'s default
# na.action = na.omit silently dropped those rows, so
# length(predict(fit, newdata)) < nrow(newdata).
#
# Fix: detect NA in factor/character newdata columns up front; error with
# a clear message naming the column(s). Mirrors the BUG-004 OOV path.

test_that("predict errors when newdata has NA in a factor column", {
  set.seed(1)
  df <- data.frame(num = stats::runif(50),
                   cat = factor(sample(c("a", "b"), 50, replace = TRUE)))
  yy <- df$num + (df$cat == "b") + stats::rnorm(50, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$cat[1:5] <- NA
  expect_error(predict(fit, df2),
               "NA value")
})

test_that("predict errors when newdata has NA in a character column", {
  set.seed(2)
  df <- data.frame(num = stats::runif(40),
                   cat = sample(c("A", "B", "C"), 40, replace = TRUE),
                   stringsAsFactors = FALSE)
  yy <- df$num + as.integer(df$cat == "B") + stats::rnorm(40, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$cat[c(2, 7, 11)] <- NA_character_
  expect_error(predict(fit, df2),
               "NA value")
})

test_that("error message names the offending column", {
  set.seed(3)
  df <- data.frame(num = stats::runif(30),
                   grp = factor(sample(c("x", "y"), 30, replace = TRUE)))
  yy <- df$num + (df$grp == "y") + stats::rnorm(30, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  df2 <- df
  df2$grp[1] <- NA
  msg <- tryCatch(predict(fit, df2),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "column grp")
})

test_that("predict still works when factor columns have no NAs (no regression)", {
  set.seed(4)
  df <- data.frame(num = stats::runif(40),
                   cat = factor(sample(c("a", "b"), 40, replace = TRUE)))
  yy <- df$num + (df$cat == "b") + stats::rnorm(40, sd = 0.2)
  fit <- ares(df, yy, nthreads = 2)
  p <- predict(fit, df)
  expect_equal(length(p), nrow(df))
  expect_true(all(!is.na(p)))
})
