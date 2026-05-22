# logreg() HC0-HC3 robust covariance parity vs sandwich::vcovHC.
# sandwich is NOT a declared dependency; guard at runtime only.

test_that("logreg HC0-HC3 match sandwich::vcovHC (unweighted)", {
  skip_if_not_installed("sandwich")
  # vs ~ wt + drat + gear is a well-separated-free, comfortably
  # convergent binary model (am ~ * with more covariates separates).
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat", "gear")])
  fit <- logreg(y ~ wt + drat + gear, data = df)
  g <- glm(y ~ wt + drat + gear, data = df, family = binomial())
  for (hc in c("HC0", "HC1", "HC2", "HC3")) {
    v_rr <- roadrunner:::.logreg_vcov(fit, hc)
    v_sw <- sandwich::vcovHC(g, type = hc)
    # The small (~1e-4) gap reflects logreg using the converged-mu Fisher
    # information for the bread, vs glm's second-to-last-iterate bread.
    expect_equal(unname(v_rr), unname(v_sw), tolerance = 1e-3,
                 info = hc)
  }
})

test_that("logreg HC0-HC3 match sandwich::vcovHC (weighted)", {
  skip_if_not_installed("sandwich")
  set.seed(11)
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat")])
  w <- runif(nrow(df), 0.5, 3)
  fit <- logreg(y ~ wt + drat, data = df, weights = w)
  # glm warns "non-integer #successes" for non-integer prior weights;
  # the IRLS fit and the HC covariance are unaffected.
  g <- suppressWarnings(
    glm(y ~ wt + drat, data = df, family = binomial(), weights = w))
  for (hc in c("HC0", "HC1", "HC2", "HC3")) {
    v_rr <- roadrunner:::.logreg_vcov(fit, hc)
    v_sw <- sandwich::vcovHC(g, type = hc)
    expect_equal(unname(v_rr), unname(v_sw), tolerance = 1e-3,
                 info = hc)
  }
})

test_that("summary.logreg robust SEs match sandwich-based SEs", {
  skip_if_not_installed("sandwich")
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat")])
  fit <- logreg(y ~ wt + drat, data = df)
  g <- glm(y ~ wt + drat, data = df, family = binomial())
  s_rr <- summary(fit, robust = "HC3")
  se_sw <- sqrt(diag(sandwich::vcovHC(g, type = "HC3")))
  expect_equal(unname(s_rr$coefficients[, "Std. Error"]),
               unname(se_sw), tolerance = 1e-3)
})

test_that("logreg classical vcov equals robust = 'none'", {
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat")])
  fit <- logreg(y ~ wt + drat, data = df)
  expect_identical(roadrunner:::.logreg_vcov(fit, "none"), fit$vcov)
})

test_that("logreg robust HC variants differ from the classical vcov", {
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat")])
  fit <- logreg(y ~ wt + drat, data = df)
  for (hc in c("HC0", "HC1", "HC2", "HC3")) {
    v_rr <- roadrunner:::.logreg_vcov(fit, hc)
    expect_false(isTRUE(all.equal(unname(v_rr), unname(fit$vcov))),
                 info = hc)
  }
})
