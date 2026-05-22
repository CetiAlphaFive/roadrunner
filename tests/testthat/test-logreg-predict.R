# logreg() predict: link / response / class + se.fit vs predict.glm.

test_that("predict.logreg link predictions match predict.glm", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  nd <- df[1:8, ]
  expect_equal(unname(predict(fit, nd, type = "link")),
               unname(predict(g, nd, type = "link")), tolerance = 1e-7)
})

test_that("predict.logreg response predictions match predict.glm", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  nd <- df[1:8, ]
  expect_equal(unname(predict(fit, nd, type = "response")),
               unname(predict(g, nd, type = "response")),
               tolerance = 1e-7)
})

test_that("predict.logreg with newdata = NULL returns fitted values", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  expect_equal(predict(fit, type = "response"), fit$fitted.values,
               tolerance = 1e-10)
  expect_equal(predict(fit, type = "link"), fit$linear.predictors,
               tolerance = 1e-10)
})

test_that("predict.logreg type = 'class' thresholds the probability", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  nd <- df[1:10, ]
  cls <- predict(fit, nd, type = "class")
  prob <- predict(fit, nd, type = "response")
  expect_equal(cls, as.integer(prob > 0.5))
  # A higher threshold yields no more positive predictions.
  cls_hi <- predict(fit, nd, type = "class", threshold = 0.9)
  expect_true(sum(cls_hi) <= sum(cls))
})

test_that("predict.logreg type = 'class' returns factor for a factor fit", {
  df <- data.frame(y = factor(ifelse(mtcars$am > 0, "auto", "manual")),
                   mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  cls <- predict(fit, df[1:6, ], type = "class")
  expect_s3_class(cls, "factor")
  expect_equal(levels(cls), levels(df$y))
})

test_that("predict.logreg se.fit (link) matches predict.glm se.fit", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  nd <- df[1:10, ]
  p_rr <- predict(fit, nd, type = "link", se.fit = TRUE)
  p_g <- predict(g, nd, type = "link", se.fit = TRUE)
  expect_equal(unname(attr(p_rr, "se.fit")), unname(p_g$se.fit),
               tolerance = 1e-3)
})

test_that("predict.logreg se.fit (response) matches predict.glm se.fit", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  nd <- df[1:10, ]
  p_rr <- predict(fit, nd, type = "response", se.fit = TRUE)
  p_g <- predict(g, nd, type = "response", se.fit = TRUE)
  expect_equal(unname(attr(p_rr, "se.fit")), unname(p_g$se.fit),
               tolerance = 1e-3)
})

test_that("predict.logreg rejects se.fit for type = 'class'", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  expect_error(predict(fit, df[1:3, ], type = "class", se.fit = TRUE),
               "se.fit")
})

test_that("predict.logreg matrix newdata works", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- as.integer(mtcars$am)
  fit <- logreg(x, y)
  p <- predict(fit, x[1:5, ], type = "response")
  expect_length(p, 5L)
  expect_equal(unname(p), unname(fit$fitted.values[1:5]),
               tolerance = 1e-10)
})

test_that("predict.logreg rejects newdata with wrong column count", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  fit <- logreg(x, as.integer(mtcars$am))
  expect_error(predict(fit, x[, 1, drop = FALSE]), "columns")
})

test_that("predict.logreg robust se.fit differs from classical", {
  skip_if_not_installed("sandwich")
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  nd <- df[1:10, ]
  se_cl <- attr(predict(fit, nd, type = "link", se.fit = TRUE), "se.fit")
  se_rb <- attr(predict(fit, nd, type = "link", se.fit = TRUE,
                        robust = "HC3"), "se.fit")
  expect_false(isTRUE(all.equal(se_cl, se_rb)))
})
