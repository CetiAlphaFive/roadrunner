# roadrunner (development)

## Package

- Renamed from `ares` to `roadrunner`. The MARS fitter remains
  available as `ares()`, with the `"ares"` S3 class and methods
  unchanged. Update `library(ares)` calls to `library(roadrunner)`.

## `ares()`

- `family` accepts `"gaussian"` (default), `"binomial"`, `"poisson"`,
  and `"gamma"`. The forward + backward MARS pass runs on the
  numeric response; the selected basis is then refit with
  `stats::glm.fit()` for non-gaussian families. `predict()` gains
  `type = c("response", "link")`.
- `weights` argument for observation weights. Composes with CV
  pruning, autotune, and bagging.
- `pmethod = "cv"` (and the convenience trigger `nfold > 0`) performs
  K-fold cross-validated subset-size selection. Optional `ncross`
  repetitions, quantile-based stratification, and the 1-SE rule via
  `cv.1se = TRUE`.
- `autotune = TRUE` runs an inner-CV grid search over
  `(degree, penalty, nk, fast.k)`. `autotune.speed` chooses among
  `"balanced"` (default), `"quality"`, and `"fast"`. A 20 % subsample
  warm-start short-circuits the full grid when the best-per-degree gap
  is decisive.
- `n.boot > 0` fits row-bootstrap replicates of the central model.
  `predict()` averages across replicates; `se.fit = TRUE` attaches the
  per-row bag standard deviation.
- `varmod = "const" | "lm"` (gaussian only) stores a residual variance
  model at fit time, enabling `predict(interval = "pint")` for
  approximate prediction intervals.
- `na.action = c("impute", "omit")` handles missing values in `x`.
  Default `"impute"` median-imputes numeric columns and stores the
  medians for reapplication at predict time.
- Factor and character columns in a data-frame `x` are expanded via
  `model.matrix` and replayed on new data.

## Engine

- Fits are parallel-deterministic: at a fixed `seed.cv`, results are
  byte-identical across thread counts.
- Default `auto.linpreds = TRUE` and `adjust.endspan = 2L` for
  `family = "binomial"` (matches `earth`'s binomial defaults). Other
  families keep the gaussian-conservative defaults.

## Compatibility

- Requires R (>= 4.1.0).
- Depends only on `Rcpp` and `RcppParallel`. `earth`, `bench`, and
  `testthat` are Suggests only.
