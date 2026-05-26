# bgam smoke test -- inst/sims/scratch/bgam_smoke.R
# Verifies:
#   1. Fit on small simulated gaussian dataset (n=120, p=5)
#   2. Output object has exactly the 25 required fields
#   3. predict() runs on new data
#   4. nthreads=1 vs nthreads=4: fitted values are byte-identical
#   5. Binomial fit basic sanity check
#   6. print() and summary() run without error

devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)

cat("=== bgam smoke test ===\n\n")

pass <- function(msg) cat("[PASS]", msg, "\n")
fail <- function(msg) { cat("[FAIL]", msg, "\n"); stop(msg) }

# ---- Test 1: Gaussian fit -----------------------------------------------
set.seed(42L)
n <- 120L
x <- matrix(rnorm(n * 5), n, 5,
             dimnames = list(NULL, paste0("x", 1:5)))
y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + 0.3 * x[, 3] + rnorm(n, sd = 0.5)

fit <- bgam(x, y, mstop = 50L, autotune = FALSE, nthreads = 1L)

if (!inherits(fit, "bgam")) fail("fit not of class bgam")
pass("Test 1a: fit has class 'bgam'")

# ---- Test 2: 25-field output contract ------------------------------------
required_fields <- c(
  "coefficients", "fitted.values", "linear.predictors", "residuals",
  "selection_path", "selection_frequency", "loss_path", "nu", "mstop",
  "mstop_opt", "family", "base_learners", "nknots", "degree", "dpen",
  "df_target", "weights", "n", "p", "predictor_names", "intercept_value",
  "sigma2", "cv", "boot", "call"
)
# terms, xlevels, na.action are also present but only guaranteed for formula
# path or when NAs exist; they should be present as NULL fields
required_fields_ext <- c(required_fields, "terms", "xlevels", "na.action")

missing_fields <- setdiff(required_fields, names(fit))
if (length(missing_fields) > 0L)
  fail(paste("Missing fields:", paste(missing_fields, collapse = ", ")))
pass("Test 2: all 25 required fields present")

# Verify specific fields
if (length(fit$fitted.values) != n) fail("fitted.values length mismatch")
if (length(fit$selection_path) != 50L) fail("selection_path length != mstop")
if (length(fit$loss_path) != 50L) fail("loss_path length != mstop")
if (fit$mstop != 50L) fail("mstop field wrong")
if (!is.null(fit$mstop_opt)) fail("mstop_opt should be NULL when autotune=FALSE")
if (fit$family != "gaussian") fail("family wrong")
if (fit$n != n) fail("n field wrong")
if (fit$p != 5L) fail("p field wrong")
if (length(fit$base_learners) != 5L) fail("base_learners length wrong")
if (length(fit$selection_frequency) != 5L) fail("selection_frequency length wrong")
if (abs(sum(fit$selection_frequency) - 1) > 1e-12)
  fail("selection_frequency doesn't sum to 1")
pass("Test 2b: output contract fields have correct types/lengths")

# ---- Test 3: predict() on new data --------------------------------------
x_new <- matrix(rnorm(10 * 5), 10, 5)
p_resp <- predict(fit, x_new, type = "response")
p_link <- predict(fit, x_new, type = "link")
p_null <- predict(fit)

if (length(p_resp) != 10L) fail("predict response length wrong")
if (length(p_link) != 10L) fail("predict link length wrong")
if (length(p_null) != n)   fail("predict NULL (training) length wrong")
if (!all(is.finite(p_resp))) fail("predict response has non-finite values")
# For gaussian, link == response
if (max(abs(p_resp - p_link)) > 1e-12) fail("gaussian link != response")
pass("Test 3: predict() on new data and NULL")

# type = "terms" returns n x p matrix
p_terms <- predict(fit, x_new, type = "terms")
if (!is.matrix(p_terms) || ncol(p_terms) != 5L || nrow(p_terms) != 10L)
  fail("predict type=terms wrong dimensions")
pass("Test 3b: predict(type='terms') returns correct matrix")

# ---- Test 4: nthreads=1 vs nthreads=4 byte-identical -------------------
fit1 <- bgam(x, y, mstop = 50L, autotune = FALSE, nthreads = 1L)
fit4 <- bgam(x, y, mstop = 50L, autotune = FALSE, nthreads = 4L)

max_diff_fitted  <- max(abs(fit1$fitted.values - fit4$fitted.values))
max_diff_selpath <- max(abs(fit1$selection_path - fit4$selection_path))

if (max_diff_fitted > 0) fail(paste0("nthreads=1 vs 4: fitted differ by ",
                                     max_diff_fitted))
if (max_diff_selpath > 0) fail("nthreads=1 vs 4: selection_path differs")
pass("Test 4: nthreads=1 vs nthreads=4 byte-identical (gaussian)")

# ---- Test 5: Binomial fit -----------------------------------------------
y_bin <- as.numeric(y > median(y))
fit_bin <- bgam(x, y_bin, family = "binomial", mstop = 30L,
                autotune = FALSE, nthreads = 1L)

if (fit_bin$family != "binomial") fail("binomial family field wrong")
if (!all(fit_bin$fitted.values >= 0 & fit_bin$fitted.values <= 1))
  fail("binomial fitted values outside [0,1]")
if (!is.null(fit_bin$sigma2)) fail("sigma2 should be NULL for binomial")

# Test nthreads=1 vs 4 for binomial
fit_bin1 <- bgam(x, y_bin, family = "binomial", mstop = 30L,
                 autotune = FALSE, nthreads = 1L)
fit_bin4 <- bgam(x, y_bin, family = "binomial", mstop = 30L,
                 autotune = FALSE, nthreads = 4L)
max_diff_bin <- max(abs(fit_bin1$fitted.values - fit_bin4$fitted.values))
if (max_diff_bin > 0) fail(paste0("nthreads=1 vs 4 (binomial): differ by ",
                                  max_diff_bin))
pass("Test 5: binomial fit + nthreads byte-identical")

# ---- Test 6: print and summary ------------------------------------------
out_print <- capture.output(print(fit))
if (length(out_print) == 0L) fail("print produced no output")

s <- summary(fit)
if (!inherits(s, "summary.bgam")) fail("summary wrong class")
out_summary <- capture.output(print(s))
if (length(out_summary) == 0L) fail("summary print produced no output")
pass("Test 6: print() and summary() run without error")

# ---- Test 7: formula interface ------------------------------------------
df <- as.data.frame(x)
df$y <- y
fit_f <- bgam(y ~ x1 + x2 + x3 + x4 + x5, data = df,
              mstop = 50L, autotune = FALSE, nthreads = 1L)
if (!inherits(fit_f, "bgam")) fail("formula fit not class bgam")
if (fit_f$p != 5L) fail("formula fit p wrong")
p_formula <- predict(fit_f, df[1:10, ])
if (length(p_formula) != 10L) fail("formula predict length wrong")
pass("Test 7: formula interface + predict")

cat("\n=== All smoke tests PASSED ===\n")
