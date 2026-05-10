test_that("ares fits and returns expected structure on mtcars", {
  fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
  expect_s3_class(fit, "ares")
  expect_true(all(c("coefficients", "bx", "dirs", "cuts", "selected.terms",
                    "rss", "gcv", "rss.per.subset", "gcv.per.subset",
                    "fitted.values", "residuals", "namesx", "call",
                    "nk", "thresh", "penalty", "minspan", "endspan",
                    "degree", "pmethod", "nthreads") %in% names(fit)))
  expect_equal(nrow(fit$bx), length(mtcars$mpg))
  expect_equal(ncol(fit$bx), length(fit$selected.terms))
  expect_equal(length(fit$coefficients), length(fit$selected.terms))
  expect_equal(nrow(fit$dirs), nrow(fit$cuts))
  expect_equal(ncol(fit$dirs), ncol(mtcars) - 1L)
  expect_equal(fit$selected.terms[1], 1L)
})

test_that("predict returns fitted.values on training data (identity)", {
  x <- as.matrix(mtcars[, -1]); y <- mtcars$mpg
  fit <- ares(x, y, nthreads = 2)
  expect_identical(predict(fit), fit$fitted.values)
  expect_equal(predict(fit, x), fit$fitted.values, tolerance = 1e-10)
})

test_that("formula method works and matches matrix interface", {
  f1 <- ares(mpg ~ ., data = mtcars, nthreads = 2)
  f2 <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-10)
  expect_equal(f1$coefficients, f2$coefficients, tolerance = 1e-10,
               ignore_attr = TRUE)
})
