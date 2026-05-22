# ols() HC0-HC3 robust covariance parity vs sandwich::vcovHC.
# sandwich is NOT a declared dependency; guard at runtime only.

# sandwich is resolved at runtime via getExportedValue() rather than `::`
# so R CMD check does not flag an unstated dependency (mirrors the
# penalizedLDA pattern used by the plda parity tests).
.sw_vcovHC <- function() getExportedValue("sandwich", "vcovHC")

test_that("ols HC0-HC3 match sandwich::vcovHC (unweighted)", {
  skip_if_not_installed("sandwich")
  vcovHC <- .sw_vcovHC()
  fit <- ols(mpg ~ wt + hp + disp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp + disp, data = mtcars)
  for (hc in c("HC0", "HC1", "HC2", "HC3")) {
    v_rr <- roadrunner:::.ols_vcov(fit, hc)
    v_sw <- vcovHC(lm0, type = hc)
    expect_equal(unname(v_rr), unname(v_sw), tolerance = 1e-8,
                 info = hc)
  }
})

test_that("ols HC0-HC3 match sandwich::vcovHC (weighted)", {
  skip_if_not_installed("sandwich")
  vcovHC <- .sw_vcovHC()
  set.seed(11)
  w <- runif(nrow(mtcars), 0.5, 3)
  fit <- ols(mpg ~ wt + hp, data = mtcars, weights = w)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars, weights = w)
  for (hc in c("HC0", "HC1", "HC2", "HC3")) {
    v_rr <- roadrunner:::.ols_vcov(fit, hc)
    v_sw <- vcovHC(lm0, type = hc)
    expect_equal(unname(v_rr), unname(v_sw), tolerance = 1e-8,
                 info = hc)
  }
})

test_that("summary.ols robust SEs match sandwich-based SEs", {
  skip_if_not_installed("sandwich")
  vcovHC <- .sw_vcovHC()
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  s_rr <- summary(fit, robust = "HC3")
  se_sw <- sqrt(diag(vcovHC(lm0, type = "HC3")))
  expect_equal(unname(s_rr$coefficients[, "Std. Error"]),
               unname(se_sw), tolerance = 1e-8)
})

test_that("ols classical vcov equals robust = 'none'", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  expect_identical(roadrunner:::.ols_vcov(fit, "none"), fit$vcov)
})
