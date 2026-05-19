# Phase Q1 parity batch (v0.0.0.9048):
#   A3 - Neffective (effective df) on summary.krls_rr.
#   A7 - subset arg on krls.default.
#   A8 - predict(..., type = "prob") for binary fits.
#   A9 - slim_krls() / unslim_krls() helpers.
#
# All tests are quick (n <= 120) and run by default (no skip_if_quick).

make_dgp <- function(n = 80L, p = 3L, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L])
  if (p >= 2L) y <- y + 0.5 * X[, 2L]^2
  if (p >= 3L) y <- y - 0.3 * X[, 3L]
  y <- y + rnorm(n, sd = 0.2)
  list(X = X, y = y)
}

make_binary_dgp <- function(n = 80L, p = 3L, seed = 2L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  eta <- 0.7 * X[, 1L] - 0.4 * X[, 2L] + 0.2 * X[, 3L]
  y <- as.numeric(stats::plogis(eta) > runif(n))
  list(X = X, y = y)
}

# ----------------------------------------------------------------------
# A3 - Neffective
# ----------------------------------------------------------------------
test_that("A3: Neffective is computed and printed on summary", {
  d <- make_dgp(n = 60L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, derivative = FALSE,
                          binary = FALSE, vcov = FALSE)
  expect_true(!is.null(fit$Neffective))
  expect_true(is.finite(fit$Neffective))
  expect_gt(fit$Neffective, 0)
  expect_lte(fit$Neffective, nrow(d$X))

  s <- summary(fit)
  expect_true(!is.null(s$Neffective))
  expect_equal(s$Neffective, fit$Neffective)

  # Output line includes the effective df value.
  out_lines <- capture.output(print(s))
  expect_true(any(grepl("Effective df", out_lines, fixed = TRUE)))
})

test_that("A3: Neffective -> 0 as lambda -> Inf, -> n as lambda -> 0", {
  d <- make_dgp(n = 60L, p = 3L)
  # Huge lambda: shrinkage -> 0 per eigen mode.
  fit_big <- roadrunner::krls(d$X, d$y, lambda = 1e12,
                              derivative = FALSE, binary = FALSE,
                              vcov = FALSE)
  expect_lt(fit_big$Neffective, 1e-6)

  # Tiny lambda: shrinkage -> 1 per eigen mode; Neff approaches the
  # number of strictly positive eigenvalues (~ n for a full-rank K).
  fit_small <- roadrunner::krls(d$X, d$y, lambda = 1e-12,
                                derivative = FALSE, binary = FALSE,
                                vcov = FALSE)
  expect_gt(fit_small$Neffective, nrow(d$X) - 1)
  expect_lte(fit_small$Neffective, nrow(d$X) + 1e-6)
})

test_that("A3: Neffective matches hand-computed sum(d/(d+lambda))", {
  d <- make_dgp(n = 60L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, derivative = FALSE,
                          binary = FALSE, vcov = FALSE)
  K <- fit$K
  expect_false(is.null(K))
  evs <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  evs[evs < 0] <- 0
  hand <- sum(evs / (evs + fit$lambda))
  expect_equal(fit$Neffective, hand, tolerance = 1e-8)
})

# ----------------------------------------------------------------------
# A7 - subset
# ----------------------------------------------------------------------
test_that("A7: integer subset matches manual slicing byte-equal", {
  d <- make_dgp(n = 80L, p = 3L, seed = 11L)
  fit_sub <- roadrunner::krls(d$X, d$y, subset = 1:50, sigma = 3,
                              lambda = 0.5,
                              derivative = FALSE, binary = FALSE,
                              vcov = FALSE)
  fit_man <- roadrunner::krls(d$X[1:50, , drop = FALSE], d$y[1:50],
                              sigma = 3, lambda = 0.5,
                              derivative = FALSE, binary = FALSE,
                              vcov = FALSE)
  expect_equal(as.numeric(fit_sub$coeffs), as.numeric(fit_man$coeffs))
  expect_equal(as.numeric(fit_sub$fitted), as.numeric(fit_man$fitted))
  expect_equal(fit_sub$lambda, fit_man$lambda)
})

test_that("A7: logical subset matches integer subset", {
  d <- make_dgp(n = 80L, p = 3L, seed = 12L)
  log_keep <- rep(FALSE, 80L)
  log_keep[1:50] <- TRUE
  fit_log <- roadrunner::krls(d$X, d$y, subset = log_keep, sigma = 3,
                              lambda = 0.5,
                              derivative = FALSE, binary = FALSE,
                              vcov = FALSE)
  fit_int <- roadrunner::krls(d$X, d$y, subset = 1:50, sigma = 3,
                              lambda = 0.5,
                              derivative = FALSE, binary = FALSE,
                              vcov = FALSE)
  expect_equal(as.numeric(fit_log$coeffs), as.numeric(fit_int$coeffs))
})

test_that("A7: subset = NULL is byte-equal to no subset", {
  d <- make_dgp(n = 60L, p = 3L, seed = 13L)
  fit_a <- roadrunner::krls(d$X, d$y, subset = NULL, sigma = 3,
                            lambda = 0.5,
                            derivative = FALSE, binary = FALSE,
                            vcov = FALSE)
  fit_b <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                            derivative = FALSE, binary = FALSE,
                            vcov = FALSE)
  expect_equal(as.numeric(fit_a$coeffs), as.numeric(fit_b$coeffs))
  expect_equal(as.numeric(fit_a$fitted), as.numeric(fit_b$fitted))
})

test_that("A7: invalid subset errors clearly", {
  d <- make_dgp(n = 30L, p = 2L)
  expect_error(roadrunner::krls(d$X, d$y, subset = c(1, 5, 999)),
               "out of bounds")
  expect_error(roadrunner::krls(d$X, d$y,
                                subset = c(TRUE, FALSE)),
               "length nrow")
  expect_error(roadrunner::krls(d$X, d$y, subset = "all"),
               "logical or integer")
})

# ----------------------------------------------------------------------
# A8 - predict(..., type = "prob") for binary fits
# ----------------------------------------------------------------------
test_that("A8: binary fit predict(type='prob') returns probabilities", {
  d <- make_binary_dgp(n = 80L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                          derivative = FALSE, binary = FALSE,
                          vcov = FALSE)
  expect_true(isTRUE(fit$binary_y))

  pr_resp <- predict(fit, newdata = d$X)
  pr_prob <- predict(fit, newdata = d$X, type = "prob")
  expect_true(all(pr_prob$fit >= 0 & pr_prob$fit <= 1))
  expect_equal(as.numeric(pr_prob$fit),
               as.numeric(stats::plogis(pr_resp$fit)))
})

test_that("A8: non-binary fit rejects type='prob'", {
  d <- make_dgp(n = 60L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                          derivative = FALSE, binary = FALSE,
                          vcov = FALSE)
  expect_false(isTRUE(fit$binary_y))
  expect_error(predict(fit, newdata = d$X, type = "prob"),
               "requires a binary fit")
})

test_that("A8: type='link' is a synonym for 'response'", {
  d <- make_binary_dgp(n = 60L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                          derivative = FALSE, binary = FALSE,
                          vcov = FALSE)
  pr_link <- predict(fit, newdata = d$X, type = "link")
  pr_resp <- predict(fit, newdata = d$X, type = "response")
  expect_equal(as.numeric(pr_link$fit), as.numeric(pr_resp$fit))
})

test_that("A8: NULL newdata also honours type='prob'", {
  d <- make_binary_dgp(n = 60L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                          derivative = FALSE, binary = FALSE,
                          vcov = FALSE)
  pr <- predict(fit, type = "prob")
  expect_true(all(pr$fit >= 0 & pr$fit <= 1))
  expect_equal(as.numeric(pr$fit),
               as.numeric(stats::plogis(as.numeric(fit$fitted))))
})

# ----------------------------------------------------------------------
# A9 - slim_krls / unslim_krls
# ----------------------------------------------------------------------
test_that("A9: slim_krls(keep_predict=TRUE) shrinks saveRDS size", {
  d <- make_dgp(n = 120L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, derivative = TRUE, binary = FALSE,
                          vcov = TRUE)
  tmp_full <- tempfile(fileext = ".rds")
  tmp_slim <- tempfile(fileext = ".rds")
  on.exit({
    if (file.exists(tmp_full)) file.remove(tmp_full)
    if (file.exists(tmp_slim)) file.remove(tmp_slim)
  }, add = TRUE)
  saveRDS(fit, tmp_full)
  saveRDS(slim_krls(fit), tmp_slim)
  sz_full <- file.info(tmp_full)$size
  sz_slim <- file.info(tmp_slim)$size
  expect_lt(sz_slim, 0.5 * sz_full)
})

test_that("A9: predict still works after slim with keep_predict=TRUE", {
  d <- make_dgp(n = 80L, p = 3L)
  fit  <- roadrunner::krls(d$X, d$y, sigma = 3, lambda = 0.5,
                           derivative = FALSE, binary = FALSE,
                           vcov = FALSE)
  pred_full <- predict(fit, newdata = d$X)$fit
  slim <- slim_krls(fit, keep_predict = TRUE)
  expect_true(isTRUE(slim$slimmed))
  expect_true(isTRUE(slim$slim_keep_predict))
  pred_slim <- predict(slim, newdata = d$X)$fit
  expect_equal(as.numeric(pred_slim), as.numeric(pred_full))
})

test_that("A9: predict errors after slim with keep_predict=FALSE", {
  d <- make_dgp(n = 80L, p = 3L)
  fit  <- roadrunner::krls(d$X, d$y, derivative = TRUE, binary = FALSE,
                           vcov = TRUE)
  slim <- slim_krls(fit, keep_predict = FALSE)
  expect_true(isTRUE(slim$slimmed))
  expect_false(isTRUE(slim$slim_keep_predict))
  expect_error(predict(slim, newdata = d$X), "slimmed")
})

test_that("A9: slim preserves Neffective (cheap field)", {
  d <- make_dgp(n = 80L, p = 3L)
  fit  <- roadrunner::krls(d$X, d$y, derivative = FALSE, binary = FALSE,
                           vcov = FALSE)
  slim_t <- slim_krls(fit, keep_predict = TRUE)
  slim_f <- slim_krls(fit, keep_predict = FALSE)
  expect_equal(slim_t$Neffective, fit$Neffective)
  expect_equal(slim_f$Neffective, fit$Neffective)
})

test_that("A9: unslim_krls is a no-op with a warning", {
  d <- make_dgp(n = 60L, p = 3L)
  fit  <- roadrunner::krls(d$X, d$y, derivative = FALSE, binary = FALSE,
                           vcov = FALSE)
  slim <- slim_krls(fit, keep_predict = TRUE)
  expect_warning(out <- unslim_krls(slim), "refit")
  expect_true(isTRUE(out$slimmed))

  # A non-slimmed fit goes through unchanged with no warning.
  expect_silent(out2 <- unslim_krls(fit))
  expect_identical(class(out2), class(fit))
})

# ----------------------------------------------------------------------
# Composition guards (per spec)
# ----------------------------------------------------------------------
test_that("A7 + A8 compose: subset on a binary DGP keeps binary detection", {
  d <- make_binary_dgp(n = 80L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, subset = 1:60, sigma = 3,
                          lambda = 0.5,
                          derivative = FALSE, binary = FALSE,
                          vcov = FALSE)
  expect_true(isTRUE(fit$binary_y))
  pr <- predict(fit, newdata = d$X, type = "prob")
  expect_true(all(pr$fit >= 0 & pr$fit <= 1))
})

test_that("A9 + A3 compose: slim keeps Neffective even with keep_predict=FALSE", {
  d <- make_dgp(n = 80L, p = 3L)
  fit <- roadrunner::krls(d$X, d$y, derivative = TRUE, binary = FALSE,
                          vcov = TRUE)
  slim <- slim_krls(fit, keep_predict = FALSE)
  expect_true(!is.null(slim$Neffective))
  expect_equal(slim$Neffective, fit$Neffective)
})
