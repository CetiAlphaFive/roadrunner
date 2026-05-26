# Simulation harness for bgam vs baselines
# REQ-20260526-013349-mboost-gam
# sim-spec.md Sections 5, 6, 9
#
# SEED DERIVATION RULE:
#   master_seed <- 20260526
#   seed_for_rep <- master_seed + dgp_index * 1000000L + cell_index * 10000L + rep_index
#   Where: dgp_index: DGP1=1, DGP2=2, DGP3=3
#          cell_index: 1-based within the DGP's cell list (see SCENARIO_GRID below)
#          rep_index: 1..200
#   RNG type: Mersenne-Twister / Inversion (set in each dgp_*_generate call)
#
# DETERMINISM NOTE: Seeds are set INSIDE dgp_generate calls, not globally.
# Each (dgp, cell, rep) has a unique deterministic seed, so replications are
# independent and results are reproducible across runs and parallelization strategies.
#
# PARALLELIZATION: mclapply over (cell, rep) pairs; per-worker seeds set inside DGP.
# On Windows: falls back to lapply (mclapply uses fork-based parallelism).

MASTER_SEED <- 20260526L

# ---------------------------------------------------------------------------
# Scenario grids (Section 4 of sim-spec)
# ---------------------------------------------------------------------------

# DGP1 grid: n x p x sigma = 2 x 2 x 2 = 8 cells
DGP1_GRID <- expand.grid(
  n     = c(200L, 1000L),
  p     = c(10L, 50L),
  sigma = c(0.5, 1.0),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

# DGP2 grid: n x p x sigma = 2 x 1 x 2 = 4 cells
DGP2_GRID <- expand.grid(
  n     = c(200L, 1000L),
  p     = c(10L),
  sigma = c(0.5, 1.0),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

# DGP3 grid: n x p = 2 x 2 = 4 cells
DGP3_GRID <- expand.grid(
  n = c(200L, 1000L),
  p = c(10L, 50L),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

N_TEST  <- 1000L
N_REPS  <- 200L

# ---------------------------------------------------------------------------
# Learner fit + predict helpers
# ---------------------------------------------------------------------------

# Build polynomial expansion for ols/logreg (stores poly object for test-set use)
build_poly_x <- function(X_train, n_signal = 3L, degree = 3L) {
  p <- ncol(X_train)
  poly_objs <- vector("list", n_signal)
  parts     <- vector("list", p)
  for (j in seq_len(n_signal)) {
    poly_objs[[j]] <- poly(X_train[, j], degree)
    parts[[j]]     <- poly_objs[[j]]
  }
  if (p > n_signal) {
    for (j in seq(n_signal + 1L, p)) {
      parts[[j]] <- X_train[, j, drop = FALSE]
    }
  }
  Xp <- do.call(cbind, parts)
  list(Xp = Xp, poly_objs = poly_objs)
}

apply_poly_x <- function(X_test, poly_objs, n_signal = 3L) {
  p     <- ncol(X_test)
  parts <- vector("list", p)
  for (j in seq_len(n_signal)) {
    parts[[j]] <- predict(poly_objs[[j]], X_test[, j])
  }
  if (p > n_signal) {
    for (j in seq(n_signal + 1L, p)) {
      parts[[j]] <- X_test[, j, drop = FALSE]
    }
  }
  do.call(cbind, parts)
}

# AUC via Wilcoxon rank-sum (no pROC dependency)
compute_auc <- function(y, p_hat) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r  <- rank(p_hat)
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

# Calibration slope: logistic regression of y on logit(p_hat)
compute_cal_slope <- function(y, p_hat) {
  p_hat_clip <- pmin(pmax(p_hat, 0.001), 0.999)
  logit_p    <- log(p_hat_clip / (1 - p_hat_clip))
  if (stats::var(logit_p) < 1e-12) return(NA_real_)
  tryCatch({
    fit_cal <- stats::glm(y ~ logit_p, family = stats::binomial())
    unname(stats::coef(fit_cal)[2L])
  }, error = function(e) NA_real_)
}

# ---------------------------------------------------------------------------
# bgam availability check
# ---------------------------------------------------------------------------
bgam_available <- function() {
  tryCatch(
    exists("bgam", envir = asNamespace("roadrunner"), inherits = FALSE),
    error = function(e) FALSE
  )
}

# Dispatch wrapper: errors informatively if bgam not yet built
dispatch_bgam <- function(X_train, y_train, family) {
  if (!bgam_available()) {
    stop(paste0(
      "bgam not yet available — rerun after builder completes. ",
      "Use run_bgam_cells.R once feat/bgam is compiled."
    ))
  }
  roadrunner::bgam(
    X_train, y_train,
    family     = family,
    autotune   = TRUE,
    mstop_max  = 200L,
    nfold      = 5L,
    seed.cv    = 0L,
    nknots     = 20L,
    degree     = 3L,
    dpen       = 2L,
    nthreads   = 1L
  )
}

# ---------------------------------------------------------------------------
# Per-rep worker: fits all learners and returns a list of metric rows
# ---------------------------------------------------------------------------
run_one_rep <- function(dgp_id, dgp_index, cell_index, cell, rep_index, learners) {

  seed_rep <- MASTER_SEED + dgp_index * 1000000L + cell_index * 10000L + rep_index

  # ---- generate data -------------------------------------------------------
  if (dgp_id == "DGP1") {
    train <- dgp1_generate(cell$n, cell$p, cell$sigma, seed_rep)
    test  <- dgp1_generate(N_TEST, cell$p, cell$sigma,
                           seed_rep + 500000L)  # distinct test seed
  } else if (dgp_id == "DGP2") {
    train <- dgp2_generate(cell$n, cell$p, cell$sigma, seed_rep)
    test  <- dgp2_generate(N_TEST, cell$p, cell$sigma,
                           seed_rep + 500000L)
  } else {
    train <- dgp3_generate(cell$n, cell$p, seed_rep)
    test  <- dgp3_generate(N_TEST, cell$p,
                           seed_rep + 500000L)
  }

  family <- if (dgp_id == "DGP3") "binomial" else "gaussian"

  # ---- common predict targets ----------------------------------------------
  X_train <- train$X
  y_train <- train$y
  X_test  <- test$X
  y_test  <- test$y

  # ---- poly expansion (for ols/logreg) -------------------------------------
  n_signal <- 3L  # poly for first 3 cols; only 2 are signal but spec says cols 1-3
  poly_built  <- build_poly_x(X_train, n_signal = n_signal, degree = 3L)
  X_poly_tr   <- poly_built$Xp
  poly_objs   <- poly_built$poly_objs
  X_poly_te   <- apply_poly_x(X_test, poly_objs, n_signal = n_signal)

  # ---- loop over learners --------------------------------------------------
  rows <- vector("list", length(learners))
  for (li in seq_along(learners)) {
    lname <- learners[li]

    base_row <- list(
      dgp      = dgp_id,
      cell_id  = cell_index,
      rep      = rep_index,
      learner  = lname,
      n        = cell$n,
      p        = cell$p,
      sigma    = if (!is.null(cell$sigma)) cell$sigma else NA_real_,
      rmse        = NA_real_,
      bias        = NA_real_,
      tp_rate     = NA_real_,
      coverage_95 = NA_real_,
      pi_width    = NA_real_,
      brier       = NA_real_,
      auc         = NA_real_,
      cal_slope   = NA_real_,
      failed      = FALSE,
      error_msg   = NA_character_
    )

    result <- tryCatch({

      if (lname == "bgam") {
        # ---- bgam --------------------------------------------------------
        fit <- dispatch_bgam(X_train, y_train, family)
        y_hat <- as.numeric(predict(fit, newdata = X_test, type = "response"))

        # TP rate: fraction of iterations selecting signal predictors (cols 1-2)
        tp_rate <- NA_real_
        if (!is.null(fit$selection_freq)) {
          sf      <- fit$selection_freq  # named vector, length p
          tp_rate <- (sf[1] + sf[2]) / sum(sf)
        }

        if (family == "gaussian") {
          rmse <- sqrt(mean((y_test - y_hat)^2))
          bias <- mean(y_hat - y_test)
          modifyList(base_row, list(rmse = rmse, bias = bias, tp_rate = tp_rate))
        } else {
          brier     <- mean((y_test - y_hat)^2)
          auc_val   <- compute_auc(y_test, y_hat)
          cal_slope <- compute_cal_slope(y_test, y_hat)
          modifyList(base_row, list(brier = brier, auc = auc_val,
                                    cal_slope = cal_slope, tp_rate = tp_rate))
        }

      } else if (lname == "ares") {
        # ---- ares --------------------------------------------------------
        if (dgp_id == "DGP2") {
          # DGP2: fit with varmod for PI
          fit   <- roadrunner::ares(X_train, y_train, family = "gaussian",
                                    autotune = TRUE, varmod = "lm",
                                    nthreads = 1L)
          pi_m  <- predict(fit, newdata = X_test, interval = "pint", level = 0.95)
          y_hat <- as.numeric(pi_m[, "fit"])
          lwr   <- as.numeric(pi_m[, "lwr"])
          upr   <- as.numeric(pi_m[, "upr"])
          cov95    <- mean(y_test >= lwr & y_test <= upr)
          pi_width <- mean(upr - lwr)
          rmse     <- sqrt(mean((y_test - y_hat)^2))
          bias     <- mean(y_hat - y_test)
          modifyList(base_row, list(rmse = rmse, bias = bias,
                                    coverage_95 = cov95, pi_width = pi_width))
        } else if (family == "gaussian") {
          fit   <- roadrunner::ares(X_train, y_train, family = "gaussian",
                                    autotune = TRUE, nthreads = 1L)
          y_hat <- as.numeric(predict(fit, newdata = X_test, type = "response"))
          rmse  <- sqrt(mean((y_test - y_hat)^2))
          bias  <- mean(y_hat - y_test)
          modifyList(base_row, list(rmse = rmse, bias = bias))
        } else {
          # binomial
          fit   <- roadrunner::ares(X_train, y_train, family = "binomial",
                                    autotune = TRUE, nthreads = 1L)
          y_hat <- as.numeric(predict(fit, newdata = X_test, type = "response"))
          brier     <- mean((y_test - y_hat)^2)
          auc_val   <- compute_auc(y_test, y_hat)
          cal_slope <- compute_cal_slope(y_test, y_hat)
          modifyList(base_row, list(brier = brier, auc = auc_val,
                                    cal_slope = cal_slope))
        }

      } else if (lname == "krls") {
        # ---- krls --------------------------------------------------------
        # ard="cheap" is the closest available ARD mode in the installed
        # package; sim-spec says ard=TRUE (full ARD, Phase 2a) which is a
        # future release. predict() returns a list; $fit is always the
        # response-scale estimate (probabilities for logistic loss).
        krls_loss <- if (family == "binomial") "logistic" else "ls"
        fit <- roadrunner::krls(
          X           = X_train,
          y           = y_train,
          whichkernel = "gaussian",
          ard         = "cheap",
          loss        = krls_loss,
          derivative  = FALSE,
          vcov        = FALSE
        )
        pred_obj <- predict(fit, newdata = X_test)
        y_hat    <- as.numeric(pred_obj$fit)

        if (family == "gaussian") {
          rmse <- sqrt(mean((y_test - y_hat)^2))
          bias <- mean(y_hat - y_test)
          modifyList(base_row, list(rmse = rmse, bias = bias))
        } else {
          brier     <- mean((y_test - y_hat)^2)
          auc_val   <- compute_auc(y_test, y_hat)
          cal_slope <- compute_cal_slope(y_test, y_hat)
          modifyList(base_row, list(brier = brier, auc = auc_val,
                                    cal_slope = cal_slope))
        }

      } else if (lname == "ols") {
        # ---- ols / logreg ------------------------------------------------
        if (family == "gaussian") {
          if (dgp_id == "DGP2") {
            fit   <- roadrunner::ols(X_poly_tr, y_train)
            pi_m  <- predict(fit, newdata = X_poly_te,
                             interval = "prediction", level = 0.95)
            y_hat    <- as.numeric(pi_m[, "fit"])
            lwr      <- as.numeric(pi_m[, "lwr"])
            upr      <- as.numeric(pi_m[, "upr"])
            cov95    <- mean(y_test >= lwr & y_test <= upr)
            pi_width <- mean(upr - lwr)
            rmse     <- sqrt(mean((y_test - y_hat)^2))
            bias     <- mean(y_hat - y_test)
            modifyList(base_row, list(rmse = rmse, bias = bias,
                                      coverage_95 = cov95, pi_width = pi_width))
          } else {
            fit   <- roadrunner::ols(X_poly_tr, y_train)
            y_hat <- as.numeric(predict(fit, newdata = X_poly_te))
            rmse  <- sqrt(mean((y_test - y_hat)^2))
            bias  <- mean(y_hat - y_test)
            modifyList(base_row, list(rmse = rmse, bias = bias))
          }
        } else {
          # binomial: use logreg
          fit   <- roadrunner::logreg(X_poly_tr, y_train)
          y_hat <- as.numeric(predict(fit, newdata = X_poly_te, type = "response"))
          brier     <- mean((y_test - y_hat)^2)
          auc_val   <- compute_auc(y_test, y_hat)
          cal_slope <- compute_cal_slope(y_test, y_hat)
          modifyList(base_row, list(brier = brier, auc = auc_val,
                                    cal_slope = cal_slope))
        }

      } else {
        stop(paste("Unknown learner:", lname))
      }

    }, error = function(e) {
      modifyList(base_row, list(failed = TRUE, error_msg = conditionMessage(e)))
    })

    rows[[li]] <- result
  }

  rows
}

# ---------------------------------------------------------------------------
# Main simulation runner
# ---------------------------------------------------------------------------
run_simulation <- function(n_reps    = N_REPS,
                           learners  = c("bgam", "ares", "krls", "ols"),
                           cores     = 1L,
                           verbose   = TRUE) {

  # Source DGPs (expected to be sourced before calling this function)

  all_grids <- list(
    DGP1 = list(dgp_id = "DGP1", dgp_index = 1L, grid = DGP1_GRID),
    DGP2 = list(dgp_id = "DGP2", dgp_index = 2L, grid = DGP2_GRID),
    DGP3 = list(dgp_id = "DGP3", dgp_index = 3L, grid = DGP3_GRID)
  )

  all_rows <- list()

  for (dgp_name in names(all_grids)) {
    spec      <- all_grids[[dgp_name]]
    dgp_id    <- spec$dgp_id
    dgp_index <- spec$dgp_index
    grid      <- spec$grid
    n_cells   <- nrow(grid)

    if (verbose) {
      cat(sprintf("\n=== %s: %d cells x %d reps x %d learners ===\n",
                  dgp_id, n_cells, n_reps, length(learners)))
    }

    for (ci in seq_len(n_cells)) {
      cell <- as.list(grid[ci, ])
      if (verbose) {
        cell_desc <- paste(names(cell), unlist(cell), sep = "=", collapse = ", ")
        cat(sprintf("  [%s cell %d/%d] %s\n", dgp_id, ci, n_cells, cell_desc))
      }

      t0 <- proc.time()[["elapsed"]]

      if (cores > 1L && .Platform$OS.type == "unix") {
        # Parallel: mclapply over reps; each worker gets its own seed via dgp_generate
        rep_results <- parallel::mclapply(
          seq_len(n_reps),
          function(r) {
            run_one_rep(dgp_id, dgp_index, ci, cell, r, learners)
          },
          mc.cores  = cores,
          mc.set.seed = FALSE  # seeds are set inside run_one_rep
        )
      } else {
        rep_results <- lapply(
          seq_len(n_reps),
          function(r) {
            run_one_rep(dgp_id, dgp_index, ci, cell, r, learners)
          }
        )
      }

      # Flatten: each rep returns a list of (n_learners) rows
      for (rr in rep_results) {
        if (!is.null(rr)) {
          for (row in rr) {
            all_rows[[length(all_rows) + 1L]] <- row
          }
        }
      }

      elapsed <- proc.time()[["elapsed"]] - t0
      if (verbose) {
        cat(sprintf("    done in %.1fs\n", elapsed))
      }
    }
  }

  # Convert to data.frame
  do.call(rbind, lapply(all_rows, as.data.frame, stringsAsFactors = FALSE))
}
