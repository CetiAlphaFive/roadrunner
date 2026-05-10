test_that("print, summary, predict do not error", {
  x <- as.matrix(mtcars[, -1]); y <- mtcars$mpg
  fit <- ares(x, y, nthreads = 2)
  expect_invisible(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.ares")
  expect_invisible(print(s))
  p <- predict(fit, x[1:5, ])
  expect_length(p, 5L)
})
