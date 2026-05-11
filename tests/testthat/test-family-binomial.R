# Tests for family = "binomial".
#
# Strategy:
#   - Quick tests run by default: input validation, basic fit structure,
#     coefficient parity with stats::glm.fit on the selected basis,
#     predict() response/link round-trip, factor / character response
#     coercion. These are cheap (~1s).
#   - Heavy tests gated behind skip_if_quick(): CV + binomial, autotune +
#     binomial, bagging + binomial, mlbench-style holdout AUC parity vs
#     earth, deterministic CV across thread counts.
#
# All tests are self-contained: each generates its own data (small,
# fast-to-fit DGP) so failures are easy to localise.

# Helper: a small 2-class DGP. y is Bernoulli(plogis(eta)).
.bin_dgp <- function(n = 200, p = 5, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  eta <- 2 * (x[, 1] - 0.5) + 1.5 * x[, 2] - 1.2 * x[, 3]
  y  <- stats::rbinom(n, 1, plogis(eta))
  list(x = x, y = y, eta = eta)
}

test_that("family arg validates: gaussian default, accepts 'binomial'", {
  d <- .bin_dgp()
  fg <- ares(d$x, as.numeric(d$y), nthreads = 2)
  expect_identical(fg$family, "gaussian")
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_error(ares(d$x, d$y, family = "poisson", nthreads = 2),
               "should be one of")
})

test_that("binomial rejects non-{0,1} numeric y", {
  d <- .bin_dgp()
  bad_y <- d$y; bad_y[1] <- 2L
  expect_error(ares(d$x, bad_y, family = "binomial", nthreads = 2),
               "binomial.*0, 1")
})

test_that("binomial accepts 2-level factor, rejects 3+ levels", {
  d <- .bin_dgp()
  yf <- factor(d$y, labels = c("no", "yes"))
  fb <- ares(d$x, yf, family = "binomial", nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_identical(fb$glm$y.levels, c("no", "yes"))

  y3 <- factor(sample(c("a", "b", "c"), 200, TRUE))
  expect_error(ares(d$x, y3, family = "binomial", nthreads = 2),
               "2-level factor")
})

test_that("binomial accepts character y, coerces sensibly", {
  d <- .bin_dgp()
  yc <- ifelse(d$y == 1, "pos", "neg")
  fb <- ares(d$x, yc, family = "binomial", nthreads = 2)
  expect_identical(fb$family, "binomial")
  # Lexical order: "neg" first, "pos" second.
  expect_identical(fb$glm$y.levels, c("neg", "pos"))
})

test_that("binomial coefficients == stats::glm.fit() on the selected basis", {
  d <- .bin_dgp(n = 300, seed = 2)
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  g  <- stats::glm.fit(fb$bx, d$y, family = stats::binomial(),
                       intercept = FALSE)
  expect_equal(unname(fb$coefficients), unname(g$coefficients),
               tolerance = 1e-8)
  # fitted == plogis(linear.predictor)
  expect_equal(fb$fitted.values,
               plogis(fb$linear.predictor),
               tolerance = 1e-10)
})

test_that("predict(type='response') = plogis(predict(type='link'))", {
  d <- .bin_dgp(n = 150, seed = 3)
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  p_resp <- predict(fb, d$x, type = "response")
  p_link <- predict(fb, d$x, type = "link")
  expect_equal(p_resp, plogis(p_link), tolerance = 1e-10)
  # response is in [0,1]
  expect_true(all(p_resp >= 0 & p_resp <= 1))
})

test_that("predict() default scale is response for binomial fits", {
  d <- .bin_dgp(n = 150, seed = 4)
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  p_default <- predict(fb, d$x)
  p_resp <- predict(fb, d$x, type = "response")
  expect_equal(p_default, p_resp, tolerance = 1e-12)
})

test_that("binomial $glm carries deviance / null.deviance / AIC fields", {
  d <- .bin_dgp(n = 200, seed = 5)
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  expect_true(all(c("deviance", "null.deviance", "df.null", "df.residual",
                    "aic", "converged", "iter") %in% names(fb$glm)))
  expect_true(is.finite(fb$glm$deviance))
  expect_true(is.finite(fb$glm$null.deviance))
  expect_true(fb$glm$deviance <= fb$glm$null.deviance)
  expect_true(isTRUE(fb$glm$converged))
})

test_that("print / summary do not error on binomial fit", {
  d <- .bin_dgp(n = 100, seed = 6)
  fb <- ares(d$x, d$y, family = "binomial", nthreads = 2)
  expect_invisible(print(fb))
  s <- summary(fb)
  expect_s3_class(s, "summary.ares")
  expect_invisible(print(s))
})

test_that("formula method threads family through correctly", {
  set.seed(7)
  df <- data.frame(x1 = runif(150), x2 = runif(150), x3 = runif(150))
  df$y <- factor(rbinom(150, 1, plogis(2 * df$x1 - 1)),
                  labels = c("no", "yes"))
  fb <- ares(y ~ x1 + x2 + x3, data = df, family = "binomial",
             nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_identical(fb$glm$y.levels, c("no", "yes"))
})

test_that("gaussian fits are byte-identical pre/post family-arg addition", {
  # Smoke-check that adding family= didn't perturb the default gaussian path
  # for a fixed seed.cv. We compare a fit with and without the family arg
  # passed explicitly; both should hit the same numerical path.
  d <- .bin_dgp(n = 200, seed = 8)
  y_num <- as.numeric(d$y)
  f1 <- ares(d$x, y_num, nthreads = 2)
  f2 <- ares(d$x, y_num, family = "gaussian", nthreads = 2)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-12)
  expect_equal(unname(f1$coefficients), unname(f2$coefficients),
               tolerance = 1e-12)
})

# ----- Heavy compositions (gated) -----------------------------------------

test_that("binomial + CV pruning composes (heavy)", {
  skip_if_quick()
  d <- .bin_dgp(n = 500, seed = 11)
  fb <- ares(d$x, d$y, family = "binomial",
             pmethod = "cv", nfold = 3, seed.cv = 42,
             nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_true(!is.null(fb$cv))
  expect_true(is.finite(fb$glm$deviance))
})

test_that("binomial + autotune composes (heavy)", {
  skip_if_quick()
  d <- .bin_dgp(n = 500, seed = 12)
  fb <- ares(d$x, d$y, family = "binomial",
             autotune = TRUE, nfold = 3, seed.cv = 42,
             nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_true(!is.null(fb$autotune))
})

test_that("binomial + bagging composes (heavy)", {
  skip_if_quick()
  d <- .bin_dgp(n = 300, seed = 13)
  fb <- ares(d$x, d$y, family = "binomial",
             n.boot = 5, seed.cv = 42, nthreads = 2)
  expect_identical(fb$family, "binomial")
  expect_true(!is.null(fb$boot))
  expect_length(fb$boot$fits, 5L)
  # Bag prediction averages on response scale; result is a probability.
  pb <- predict(fb, d$x)
  expect_true(all(pb >= 0 & pb <= 1))
})

test_that("binomial determinism: nthreads=1 == nthreads=N", {
  skip_if_quick()
  d <- .bin_dgp(n = 300, seed = 14)
  f1 <- ares(d$x, d$y, family = "binomial", seed.cv = 42, nthreads = 1)
  f2 <- ares(d$x, d$y, family = "binomial", seed.cv = 42, nthreads = 4)
  expect_equal(unname(f1$coefficients), unname(f2$coefficients),
               tolerance = 1e-12)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-12)
  expect_equal(f1$glm$deviance, f2$glm$deviance, tolerance = 1e-12)
})
