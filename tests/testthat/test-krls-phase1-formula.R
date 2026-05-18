# Phase 1 tests: formula method, factor expansion, na.action, subset.

skip_if_no_krls <- function() {
  testthat::skip_if_not_installed("KRLS")
}

make_dgp_df <- function(n = 80L, seed = 1L,
                       p_num = 3L, p_cat = 1L, p_levels = 3L) {
  set.seed(seed)
  out <- list()
  for (j in seq_len(p_num)) out[[paste0("x", j)]] <- rnorm(n)
  for (j in seq_len(p_cat))
    out[[paste0("g", j)]] <- factor(
      sample(LETTERS[seq_len(p_levels)], n, replace = TRUE))
  df <- as.data.frame(out)
  # build y from whatever numeric cols exist plus first cat
  ynum <- rowSums(as.matrix(df[, paste0("x", seq_len(p_num)), drop = FALSE]))
  ycat <- if (p_cat > 0L) as.numeric(df$g1) else 0
  df$y <- ynum + ycat + rnorm(n, sd = 0.2)
  return(df)
}

# --- 1. matrix interface back-compat (byte identical) -----------------
test_that("Phase 1: matrix interface stays byte-identical pre/post", {
  set.seed(11L)
  X <- matrix(rnorm(60 * 3), 60, 3)
  colnames(X) <- c("a", "b", "c")
  y <- sin(X[, 1]) + 0.5 * X[, 2] + rnorm(60, sd = 0.1)
  fit <- krls(X, y, sigma = 3, lambda = 0.5,
              derivative = FALSE, vcov = TRUE)
  # smoke: object class, key fields
  expect_s3_class(fit, "krls_rr")
  expect_true(inherits(fit, "krls"))
  expect_equal(length(fit$coeffs), 60L)
  expect_equal(dim(fit$K), c(60L, 60L))
  # Verify the call is captured (Phase 8).
  expect_true(!is.null(fit$call))
})

# --- 2. formula interface ---------------------------------------------
test_that("Phase 1: formula method produces equivalent fit to matrix", {
  set.seed(12L)
  df <- make_dgp_df(n = 60L, p_num = 3L, p_cat = 0L)
  # matrix interface
  X <- as.matrix(df[, c("x1", "x2", "x3")])
  fm <- krls(y ~ x1 + x2 + x3, data = df,
             sigma = 3, lambda = 0.5,
             derivative = FALSE, vcov = TRUE)
  fx <- krls(X, df$y,
             sigma = 3, lambda = 0.5,
             derivative = FALSE, vcov = TRUE)
  expect_equal(as.numeric(fm$coeffs), as.numeric(fx$coeffs),
               tolerance = 1e-12)
  expect_equal(as.numeric(fm$fitted), as.numeric(fx$fitted),
               tolerance = 1e-12)
  expect_equal(fm$lambda, fx$lambda)
  # call / terms / xlevels available on formula fit
  expect_true(!is.null(fm$call))
  expect_true(!is.null(fm$terms))
})

# --- 3. factor expansion (formula path) -------------------------------
test_that("Phase 1: formula factor expansion matches manual model.matrix", {
  set.seed(13L)
  df <- make_dgp_df(n = 60L, p_num = 2L, p_cat = 1L)
  # manual expansion
  mm <- model.matrix(~ x1 + x2 + g1, data = df)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  fm_mat <- krls(mm, df$y, sigma = ncol(mm), lambda = 0.5,
                 derivative = FALSE, vcov = TRUE)
  fm_form <- krls(y ~ x1 + x2 + g1, data = df,
                  sigma = ncol(mm), lambda = 0.5,
                  derivative = FALSE, vcov = TRUE)
  expect_equal(as.numeric(fm_form$coeffs), as.numeric(fm_mat$coeffs),
               tolerance = 1e-10)
  expect_equal(as.numeric(fm_form$fitted), as.numeric(fm_mat$fitted),
               tolerance = 1e-10)
  expect_equal(colnames(fm_form$X), colnames(mm))
})

# --- 4. data.frame interface factor expansion (via krls.default) ------
test_that("Phase 1: data.frame interface expands factors", {
  set.seed(14L)
  df <- make_dgp_df(n = 60L, p_num = 2L, p_cat = 1L)
  X_df <- df[, c("x1", "x2", "g1")]
  fit <- krls(X_df, df$y, sigma = 4, lambda = 0.4,
              derivative = FALSE, vcov = TRUE)
  expect_true(!is.null(fit$factor_info))
  expect_equal(fit$factor_info$xlevels$g1, levels(df$g1))
})

# --- 5. subset honoured ----------------------------------------------
test_that("Phase 1: subset arg restricts rows for formula method", {
  set.seed(15L)
  df <- make_dgp_df(n = 80L, p_num = 3L, p_cat = 0L)
  fit_full <- krls(y ~ x1 + x2 + x3, data = df,
                   sigma = 3, lambda = 0.4,
                   derivative = FALSE, vcov = TRUE)
  fit_sub  <- krls(y ~ x1 + x2 + x3, data = df, subset = 1:40,
                   sigma = 3, lambda = 0.4,
                   derivative = FALSE, vcov = TRUE)
  expect_equal(nrow(fit_full$X), 80L)
  expect_equal(nrow(fit_sub$X), 40L)
  # subset fit shouldn't equal full fit
  expect_false(isTRUE(all.equal(as.numeric(fit_full$coeffs[1:5]),
                                as.numeric(fit_sub$coeffs[1:5]))))
})

# --- 6. na.action = "impute" -----------------------------------------
test_that("Phase 1: na.action = 'impute' median-fills X NAs", {
  set.seed(16L)
  X <- matrix(rnorm(50 * 3), 50, 3)
  colnames(X) <- c("a", "b", "c")
  y <- rowSums(X) + rnorm(50, sd = 0.1)
  X[5L, 1L] <- NA
  X[10L, 2L] <- NA
  # auto-impute (default) emits a warning
  expect_warning(
    fit <- krls(X, y, sigma = 3, lambda = 0.3,
                derivative = FALSE, vcov = TRUE),
    "median-imputed"
  )
  # n preserved, medians stored
  expect_equal(nrow(fit$X), 50L)
  expect_true(!is.null(fit$na.medians))
  # imputed values match manual median fill
  X_man <- X
  X_man[5L, 1L] <- median(X[, 1L], na.rm = TRUE)
  X_man[10L, 2L] <- median(X[, 2L], na.rm = TRUE)
  fit_man <- krls(X_man, y, sigma = 3, lambda = 0.3,
                  derivative = FALSE, vcov = TRUE)
  expect_equal(as.numeric(fit$coeffs), as.numeric(fit_man$coeffs),
               tolerance = 1e-12)
})

# --- 7. na.action = "omit" -------------------------------------------
test_that("Phase 1: na.action = 'omit' drops incomplete rows", {
  set.seed(17L)
  X <- matrix(rnorm(50 * 3), 50, 3)
  colnames(X) <- c("a", "b", "c")
  y <- rowSums(X) + rnorm(50, sd = 0.1)
  X[5L, 1L] <- NA
  expect_warning(
    fit <- krls(X, y, sigma = 3, lambda = 0.3,
                derivative = FALSE, vcov = TRUE,
                na.action = "omit"),
    "dropped"
  )
  expect_equal(nrow(fit$X), 49L)
})

# --- 8. OOV factor error (data.frame path) ---------------------------
test_that("Phase 1: predict on OOV factor errors clearly", {
  set.seed(18L)
  df <- make_dgp_df(n = 60L, p_num = 2L, p_cat = 1L)
  X_df <- df[, c("x1", "x2", "g1")]
  fit <- krls(X_df, df$y, sigma = 3, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  bad <- df[1:5, c("x1", "x2", "g1")]
  bad$g1 <- as.character(bad$g1)
  bad$g1[1] <- "ZZZ_OOV"
  expect_error(predict(fit, newdata = bad),
               "out-of-vocabulary")
})

# --- 9. OOV factor error (formula path) ------------------------------
test_that("Phase 1: predict on OOV factor errors clearly via formula", {
  set.seed(19L)
  df <- make_dgp_df(n = 60L, p_num = 2L, p_cat = 1L)
  df$g1 <- as.character(df$g1)
  fit <- krls(y ~ x1 + x2 + g1, data = df,
              sigma = 3, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  bad <- df[1:5, ]
  bad$g1[1] <- "ZZZ_OOV"
  expect_error(predict(fit, newdata = bad),
               "new level|OOV|out-of-vocabulary")
})

# --- 10. predict re-imputes NA in matrix newdata ---------------------
test_that("Phase 1: predict reapplies stored medians to NA in newdata", {
  set.seed(20L)
  X <- matrix(rnorm(60 * 2), 60, 2)
  colnames(X) <- c("a", "b")
  y <- X[, 1] + 0.3 * X[, 2] + rnorm(60, sd = 0.1)
  X[1L, 1L] <- NA
  expect_warning(
    fit <- krls(X, y, sigma = 2, lambda = 0.3,
                derivative = FALSE, vcov = TRUE),
    "median-imputed"
  )
  Xnew <- X[1:5, , drop = FALSE]
  Xnew[2, 2] <- NA
  expect_warning(pr <- predict(fit, newdata = Xnew),
                 "median-imputed")
  expect_equal(nrow(pr$fit), 5L)
})

# --- 11. predict NULL newdata returns fitted -------------------------
test_that("Phase 1: predict(fit, NULL) returns the training fitted", {
  set.seed(21L)
  X <- matrix(rnorm(40 * 2), 40, 2)
  y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  pr <- predict(fit, newdata = NULL)
  expect_equal(as.numeric(pr$fit), as.numeric(fit$fitted))
})
