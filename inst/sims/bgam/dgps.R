# DGPs for bgam vs baselines simulation
# REQ-20260526-013349-mboost-gam
# sim-spec.md Section 3
#
# Three DGPs:
#   DGP1 — Sparse Smooth Additive (Continuous)
#   DGP2 — Heteroskedastic Noise (Continuous)
#   DGP3 — Binary Outcome with Smooth Log-Odds
#
# All DGPs accept (n, params, seed) and return list(X, y, y_true, f_true).
# Seed is set at the TOP of each call (pure function: same inputs -> same output).

# ---------------------------------------------------------------------------
# DGP1 — Sparse Smooth Additive
#   Y = sin(2*pi*X1) + 0.5*X2^2 + epsilon
#   X_j ~ U(-2, 2), j = 1..p
#   epsilon ~ N(0, sigma^2)
#   Signal predictors: j in {1, 2}
# ---------------------------------------------------------------------------
dgp1_generate <- function(n, p, sigma, seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  X <- matrix(stats::runif(n * p, -2, 2), nrow = n, ncol = p)
  f_true <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2
  eps <- stats::rnorm(n, 0, sigma)
  y <- f_true + eps
  list(X = X, y = y, f_true = f_true, n = n, p = p, sigma = sigma)
}

# ---------------------------------------------------------------------------
# DGP2 — Heteroskedastic Noise
#   Y = sin(2*pi*X1) + 0.5*X2^2 + epsilon(X1)
#   epsilon(X1) ~ N(0, (1 + |X1|)^2 * sigma^2)
#   X_j ~ U(-2, 2), j = 1..p
#   Signal predictors: j in {1, 2}
# ---------------------------------------------------------------------------
dgp2_generate <- function(n, p, sigma, seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  X <- matrix(stats::runif(n * p, -2, 2), nrow = n, ncol = p)
  f_true <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2
  het_sd <- (1 + abs(X[, 1])) * sigma
  eps <- stats::rnorm(n, 0, het_sd)
  y <- f_true + eps
  list(X = X, y = y, f_true = f_true, n = n, p = p, sigma = sigma)
}

# ---------------------------------------------------------------------------
# DGP3 — Binary Outcome with Smooth Log-Odds
#   logit(p(X)) = sin(pi*X1) + X2 - 1
#   X_j ~ N(0, 1), j = 1..p
#   Signal predictors: j in {1, 2}
# ---------------------------------------------------------------------------
dgp3_generate <- function(n, p, seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  X <- matrix(stats::rnorm(n * p, 0, 1), nrow = n, ncol = p)
  eta_true <- sin(pi * X[, 1]) + X[, 2] - 1
  p_true <- 1 / (1 + exp(-eta_true))
  y <- stats::rbinom(n, 1, p_true)
  list(X = X, y = y, p_true = p_true, eta_true = eta_true,
       n = n, p = p)
}
