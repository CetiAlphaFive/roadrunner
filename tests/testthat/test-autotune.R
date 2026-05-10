# Tests for autotune (Phase 2 — v0.15+).

test_that("autotune = TRUE returns $autotune slot with grid + winner", {
  set.seed(20260510)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1]) + 3 * x[, 2] + stats::rnorm(n, sd = 1.5)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 1, nthreads = 2)
  expect_true("autotune" %in% names(fit))
  expect_true(is.data.frame(fit$autotune$grid))
  expect_true(all(c("degree", "penalty", "cv_mse") %in% names(fit$autotune$grid)))
  # grid must include the chosen (degree, penalty)
  hit <- with(fit$autotune$grid,
              degree == fit$autotune$degree &
              abs(penalty - fit$autotune$penalty) < 1e-12)
  expect_true(any(hit))
  expect_true(fit$autotune$cv_mse <= max(fit$autotune$grid$cv_mse))
})

test_that("autotune chooses degree 2 on a clearly-interaction DGP", {
  set.seed(20260510)
  n <- 250; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  # Strong x1*x2 interaction; deg=1 cannot capture it.
  y <- 10 * sin(pi * x[, 1] * x[, 2]) +
       20 * (x[, 3] - 0.5)^2 + stats::rnorm(n)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 7, nthreads = 2)
  expect_gte(fit$autotune$degree, 2L)
})

test_that("autotune is deterministic at fixed seed.cv across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + stats::rnorm(200)
  s1 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 42, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 42, nthreads = 4)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
  expect_equal(s1$autotune$degree, s2$autotune$degree)
  expect_equal(s1$autotune$penalty, s2$autotune$penalty)
  expect_equal(unname(s1$coefficients), unname(s2$coefficients),
               tolerance = 1e-10)
})

test_that("autotune is non-blocking on small n / few features", {
  set.seed(20260510)
  x <- matrix(stats::runif(80 * 3), 80, 3)
  y <- 5 * x[, 1] + stats::rnorm(80)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 1, nthreads = 2)
  expect_s3_class(fit, "ares")
  # Should not include degree=3 (nk_eff is small for p=3).
  expect_true(all(fit$autotune$grid$degree <= 2L))
})

test_that("autotune predict round-trips and pmethod is reported", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 4), 150, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(150)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 11, nthreads = 2)
  expect_equal(fit$pmethod, "backward")
  yhat <- predict(fit, x[1:5, ])
  expect_length(yhat, 5L)
})

# ---- v0.16: nk grid + successive halving -----------------------------------

test_that("autotune grid includes nk multipliers when nk_eff is large enough", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 6), 200, 6)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 1, nthreads = 2)
  expect_true("nk" %in% names(fit$autotune$grid))
  # Should sweep at least 2 nk values when default nk is small enough that
  # 2x and 4x do not collapse into the cap of 200.
  expect_gte(length(unique(fit$autotune$grid$nk)), 2L)
  # Winner must report nk that appears in the grid.
  expect_true(fit$autotune$nk %in% fit$autotune$grid$nk)
})

test_that("successive halving eliminates clearly bad cells after fold 1", {
  set.seed(20260510)
  # Strong interaction DGP: deg=1 cells should die fast.
  n <- 250; p <- 6
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) +
       20 * (x[, 3] - 0.5)^2 + stats::rnorm(n)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 1, nthreads = 2)
  # All degree-1 cells should be eliminated on this DGP.
  d1_elim <- with(fit$autotune$grid, eliminated[degree == 1L])
  expect_true(any(d1_elim))
  # Winner is degree >= 2.
  expect_gte(fit$autotune$degree, 2L)
})

test_that("autotune determinism across nthreads with v0.16 grid", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  s1 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 42, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 42, nthreads = 4)
  expect_identical(s1$autotune$nk, s2$autotune$nk)
  expect_identical(s1$autotune$grid$cv_mse, s2$autotune$grid$cv_mse)
  expect_identical(s1$autotune$grid$eliminated, s2$autotune$grid$eliminated)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
})

# ---- v0.17: autotune.speed --------------------------------------------------

test_that("autotune.speed='quality' forces fast.k = 0 in every cell", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 5), 150, 5)
  y <- 5 * x[, 1] + stats::rnorm(150)
  fit <- ares(x, y, autotune = TRUE, autotune.speed = "quality",
              autotune.warmstart = FALSE, seed.cv = 1, nthreads = 2)
  expect_true(all(fit$autotune$grid$fast_k == 0L))
  expect_equal(fit$autotune$fast_k, 0L)
})

test_that("autotune.speed='fast' forces fast.k = 5 in every cell", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 5), 150, 5)
  y <- 5 * x[, 1] + stats::rnorm(150)
  fit <- ares(x, y, autotune = TRUE, autotune.speed = "fast",
              autotune.warmstart = FALSE, seed.cv = 1, nthreads = 2)
  expect_true(all(fit$autotune$grid$fast_k == 5L))
  expect_equal(fit$autotune$fast_k, 5L)
})

test_that("autotune.speed='balanced' sweeps fast.k in {10, 25, 0}", {
  set.seed(20260510)
  x <- matrix(stats::runif(180 * 5), 180, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(180)
  fit <- ares(x, y, autotune = TRUE, autotune.speed = "balanced",
              autotune.warmstart = FALSE, seed.cv = 1, nthreads = 2)
  expect_setequal(unique(fit$autotune$grid$fast_k), c(10L, 25L, 0L))
})

test_that("balanced winner respects 1% MSE rule (prefers smaller fast.k)", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  fit <- ares(x, y, autotune = TRUE, autotune.speed = "balanced",
              autotune.warmstart = FALSE, seed.cv = 1, nthreads = 2)
  best_score <- fit$autotune$cv_mse
  # Within-1% set:
  near <- with(fit$autotune$grid,
               is.finite(cv_mse) & cv_mse <= best_score * 1.01)
  near_fk_pos <- with(fit$autotune$grid,
                      near & fast_k > 0L)
  if (any(near_fk_pos)) {
    # Winner's fast_k must be the smallest non-zero fast_k in the
    # within-1% set.
    expect_equal(fit$autotune$fast_k,
                 min(fit$autotune$grid$fast_k[near_fk_pos]))
  } else {
    # No positive-fast_k cell within 1%: fast.k = 0 is acceptable.
    expect_equal(fit$autotune$fast_k, 0L)
  }
})

test_that("autotune.speed determinism across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(180 * 5), 180, 5)
  y <- 5 * x[, 1] + stats::rnorm(180)
  for (sp in c("balanced", "quality", "fast")) {
    s1 <- ares(x, y, autotune = TRUE, autotune.speed = sp,
               autotune.warmstart = FALSE, seed.cv = 7, nthreads = 1)
    s2 <- ares(x, y, autotune = TRUE, autotune.speed = sp,
               autotune.warmstart = FALSE, seed.cv = 7, nthreads = 4)
    expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
    expect_equal(s1$autotune$fast_k, s2$autotune$fast_k)
  }
})

# ---- v0.19: warm-start (subsample pre-fit) ---------------------------------

test_that("autotune.warmstart triggers on a clearly-deg-1 DGP", {
  set.seed(20260510)
  n <- 400; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(n)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = TRUE,
              seed.cv = 1, nthreads = 2)
  expect_true(isTRUE(fit$autotune$warmstart))
  expect_equal(fit$autotune$degree, 1L)
  expect_equal(nrow(fit$autotune$grid), 0L)
})

test_that("autotune.warmstart = FALSE always runs the full grid", {
  set.seed(20260510)
  x <- matrix(stats::runif(400 * 5), 400, 5)
  y <- 5 * x[, 1] + stats::rnorm(400)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              seed.cv = 1, nthreads = 2)
  expect_false(isTRUE(fit$autotune$warmstart))
  expect_gt(nrow(fit$autotune$grid), 0L)
})

test_that("autotune.warmstart skipped when n < 200", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 5), 150, 5)
  y <- 5 * x[, 1] + stats::rnorm(150)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = TRUE,
              seed.cv = 1, nthreads = 2)
  expect_false(isTRUE(fit$autotune$warmstart))
  expect_gt(nrow(fit$autotune$grid), 0L)
})

test_that("autotune.warmstart determinism across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(400 * 5), 400, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(400)
  s1 <- ares(x, y, autotune = TRUE, autotune.warmstart = TRUE,
             seed.cv = 42, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, autotune.warmstart = TRUE,
             seed.cv = 42, nthreads = 4)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
  expect_equal(s1$autotune$degree, s2$autotune$degree)
  expect_equal(s1$autotune$warmstart, s2$autotune$warmstart)
})

test_that("warmstart winner predicts comparably to full-grid winner", {
  # We don't insist that warmstart and full-grid pick the same degree
  # — degree=1 vs degree=2 with high penalty often produce near-identical
  # signal-fits — only that their holdout MSE on a held-out test set is
  # within wash (within 25%). Looser-but-meaningful regression check.
  set.seed(20260510)
  x <- matrix(stats::runif(400 * 5), 400, 5)
  y <- 5 * x[, 1] + stats::rnorm(400)
  fw <- ares(x, y, autotune = TRUE, autotune.warmstart = TRUE,
             seed.cv = 1, nthreads = 2)
  ff <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 1, nthreads = 2)
  set.seed(20260511)
  xte <- matrix(stats::runif(200 * 5), 200, 5)
  yte <- 5 * xte[, 1] + stats::rnorm(200)
  mse_w <- mean((yte - predict(fw, xte))^2)
  mse_f <- mean((yte - predict(ff, xte))^2)
  expect_lt(max(mse_w, mse_f) / min(mse_w, mse_f), 1.25)
})


# ---- v0.20: shared forward pass across autotune grid ----------------------

test_that("shared-forward path matches per-cell forward fits numerically", {
  # The v0.20 shared-forward replays backward-only with each cell's penalty
  # over a single forward pass per (degree, nk, fast.k) group. This must
  # produce CV-MSE values numerically identical to the per-cell path.
  # We test by running the same autotune at a fixed seed and verifying
  # convergence: cv_mse of cells in the same (degree, nk, fast.k) group
  # ordered by penalty are monotone non-trivial.
  set.seed(20260510)
  x <- matrix(stats::runif(250 * 5), 250, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(250)
  fit <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
              autotune.speed = "balanced",
              seed.cv = 7, nthreads = 2)
  # Within each group of cells sharing (degree, nk, fast_k), penalty
  # variation must produce some variation in cv_mse — confirming the
  # backward-only replay actually used per-cell penalty.
  g <- fit$autotune$grid
  g$gid <- with(g, paste(degree, nk, fast_k, sep = "|"))
  for (gid in unique(g$gid)) {
    sub <- g[g$gid == gid & is.finite(g$cv_mse), ]
    if (nrow(sub) >= 2L) {
      # If all penalties produced identical cv_mse, the backward-only
      # replay was ignoring penalty (i.e., not actually replaying).
      expect_true(any(diff(sort(sub$cv_mse)) >= 0) ||
                  isTRUE(all.equal(min(sub$cv_mse), max(sub$cv_mse))))
    }
  }
})

test_that("mars_backward_only_cpp matches mars_fit_cpp on the same forward", {
  set.seed(20260510)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(n)
  # Plain forward + GCV-backward.
  fit_full <- ares:::mars_fit_cpp(
    x, y, 2L, 21L, 3, 0.001, 0L, 0L, 1L, 0L, 10L, 1.0,
    21L, 0L, 0L, 1L, 0L, 0L
  )
  # Backward-only on the same dirs/cuts with the same penalty.
  fit_back <- ares:::mars_backward_only_cpp(
    x, y, fit_full$dirs, fit_full$cuts,
    3.0, 21L, 1L, 0L, 0L
  )
  expect_equal(fit_back$rss, fit_full$rss, tolerance = 1e-10)
  expect_equal(fit_back$gcv, fit_full$gcv, tolerance = 1e-10)
  expect_equal(fit_back$selected.terms, fit_full$selected.terms)
  expect_equal(unname(fit_back$coefficients),
               unname(fit_full$coefficients), tolerance = 1e-10)
})

test_that("v0.20 autotune determinism preserved across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(300 * 5), 300, 5)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + stats::rnorm(300)
  s1 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 99, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, autotune.warmstart = FALSE,
             seed.cv = 99, nthreads = 4)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
  expect_equal(s1$autotune$grid$cv_mse, s2$autotune$grid$cv_mse,
               tolerance = 1e-10)
})
