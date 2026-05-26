# Component-wise P-spline gradient boosting

Fits a smooth, strictly additive model via component-wise functional
gradient boosting with penalised B-spline (P-spline) base-learners, one
per predictor column. At each of `mstop` boosting iterations the
algorithm selects the single base-learner that most reduces the current
loss (squared error for gaussian, binomial deviance for binomial),
updates only that component by a shrinkage-scaled Cholesky-backed ridge
solve, and repeats. The result is a smooth GAM with implicit variable
selection: predictors that are never selected retain a zero coefficient
vector, and `$selection_frequency` (fraction of iterations each
predictor was chosen) serves as a natural importance metric.

## Usage

``` r
bgam(x, ...)

# S3 method for class 'formula'
bgam(
  x,
  data = NULL,
  subset = NULL,
  na.action = stats::na.omit,
  weights = NULL,
  ...,
  y = NULL
)

# Default S3 method
bgam(
  x,
  y,
  family = c("gaussian", "binomial"),
  nu = 0.1,
  mstop = NULL,
  autotune = TRUE,
  mstop_max = 300L,
  nfold = 5L,
  seed.cv = NULL,
  nknots = 20L,
  degree = 3L,
  dpen = 2L,
  df_target = 4,
  lambda_method = c("df", "fixed"),
  lambda_fixed = NULL,
  unpenalized = NULL,
  intercept = TRUE,
  weights = NULL,
  n.boot = 0L,
  seed = NULL,
  nthreads = 1L,
  ...
)
```

## Arguments

- x:

  A numeric matrix or data frame of predictors (default method), or a
  model formula (formula method).

- ...:

  Passed to the default method from the formula method; otherwise
  currently unused.

- data:

  A data frame (formula method only).

- subset:

  Optional integer or logical row-subsetting vector passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  (formula method only).

- na.action:

  NA-handling function passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html).
  Default [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html).
  The default method drops rows with any `NA` in `x` and emits a
  warning; rows with `NA` in `y` are always rejected.

- weights:

  Optional strictly positive numeric weight vector of length `nrow(x)`.
  For `"gaussian"`, weights scale the pseudo-residual sum of squares
  used for base-learner selection and the Cholesky prefactor. For
  `"binomial"`, weights scale the IRLS working weights
  \\w_i^{(\text{eff})} = w_i \mu_i (1 - \mu_i)\\. Default `NULL`
  (unweighted).

- y:

  Response vector of length `nrow(x)`. For `family = "binomial"`, must
  be binary: 0/1 numeric, `logical`, or a two-level `factor` (second
  level = 1).

- family:

  Response family: `"gaussian"` (default, squared-error loss) or
  `"binomial"` (logistic / binomial deviance, binary response).

- nu:

  Boosting shrinkage / learning rate. Must be in `(0, 1]`. Smaller
  values require more iterations but produce smoother solution paths.
  Default `0.1`.

- mstop:

  Number of boosting iterations to run. If `NULL` and
  `autotune = FALSE`, defaults to `mstop_max`. If `autotune = TRUE`,
  `mstop` is set to the CV-selected `mstop_opt`.

- autotune:

  If `TRUE` (default), selects `mstop` by k-fold cross-validation over
  the grid `seq(10, mstop_max, by = 10)` using MSE (gaussian) or
  binomial deviance (binomial) as the criterion. A successive-halving
  early-stop rule terminates the grid scan once CV loss rises more than
  2% of its range above the running best.

- mstop_max:

  Maximum iterations for the CV grid (and the default when
  `autotune = FALSE` and `mstop = NULL`). Reduced automatically to
  `max(50, 2 * n)` when `n < 100`. Default `300L`.

- nfold:

  Number of CV folds. Default `5L`.

- seed.cv:

  Integer seed for CV fold construction. Default `NULL`. Binomial CV
  folds are stratified by response class.

- nknots:

  Number of interior B-spline knots per predictor. Knots are placed at
  equally-spaced quantiles of the predictor column. Default `20L`,
  giving `K_j = nknots + degree + 1 = 24` columns per base-learner at
  defaults.

- degree:

  B-spline polynomial degree (1 to 5). Default `3L` (cubic B-splines).

- dpen:

  Order of the difference penalty (1 to `degree`). Default `2L`
  (second-order penalty, analogous to a thin-plate-spline roughness
  penalty; Eilers & Marx 1996).

- df_target:

  Target effective degrees of freedom per base-learner, i.e.
  `tr(S_j(lambda_j))`. Must be `> 1`. Default `4`. The ridge penalty
  `lambda_j` is calibrated per predictor via
  [`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html).

- lambda_method:

  `"df"` (default): calibrate `lambda_j` from `df_target` per predictor.
  `"fixed"`: supply `lambda_j` directly via `lambda_fixed`.

- lambda_fixed:

  Numeric scalar or length-p vector of lambda values. Required when
  `lambda_method = "fixed"`; ignored otherwise.

- unpenalized:

  Optional character vector of predictor column names to enter as linear
  (unpenalised, no-spline) base-learners. Useful for known-linear
  covariates or binary dummies.

- intercept:

  If `TRUE` (default), initialise the additive predictor at `mean(y)`
  (gaussian) or `logit(clamp(mean(y), 0.001, 0.999))` (binomial). If
  `FALSE`, initialise at 0.

- n.boot:

  Number of bootstrap replicates for bagging. Default `0L` (no bagging).
  When `> 0`, [`predict()`](https://rdrr.io/r/stats/predict.html)
  returns the bag mean across the central fit and surviving replicates.
  Hyperparameters are frozen at the central-fit values; no re-autotune
  per replicate.

- seed:

  Optional integer RNG seed for the bagging loop. Default `NULL`.

- nthreads:

  TBB thread count for the parallel predictor scan per boosting
  iteration. `0` = use all available TBB threads. `1L` (default) =
  serial execution. Results are byte-identical across all `nthreads`
  values.

## Value

An S3 object of class `"bgam"`, a named list with the following key
fields (full contract: 25+ fields):

- `coefficients`:

  List of p accumulated (nu-scaled) beta vectors, one per predictor.
  Zero-vector if the predictor was never selected.

- `fitted.values`:

  Fitted values on the response scale (probabilities for binomial;
  identical to `linear.predictors` for gaussian). For bagged fits, this
  is the bag mean.

- `linear.predictors`:

  Fitted values on the link scale.

- `residuals`:

  Residuals on the response scale (gaussian: `y - fitted`; binomial:
  deviance residuals).

- `selection_path`:

  Integer vector of length `mstop`, 1-based predictor index selected at
  each iteration.

- `selection_frequency`:

  Named numeric vector of length p; fraction of iterations each
  predictor was selected. Sums to 1.

- `loss_path`:

  Numeric vector of length `mstop`; training loss (MSE or deviance)
  after each iteration.

- `mstop`:

  Integer; iterations run (equals `mstop_opt` when `autotune = TRUE`).

- `cv`:

  List with `$grid`, `$cv_loss`, `$mstop_opt` when `autotune = TRUE`;
  `NULL` otherwise.

- `boot`:

  List with `$n.boot` and `$fits` when `n.boot > 0`; `NULL` otherwise.

- `base_learners`:

  List of p per-predictor lists, each containing `$knots`, `$K`,
  `$degree`, `$dpen`, `$lambda`, `$range`, `$chol`, `$B_train`.

See also: `predict.bgam`, `print.bgam`, `summary.bgam`, `plot.bgam`.

## Details

Each base-learner is a cubic (or general degree-`d`) B-spline with
`nknots` interior knots placed at equally-spaced quantiles of the
predictor, penalised by a `dpen`-th order difference-penalty matrix
\\D^T D\\ scaled by a ridge parameter \\\lambda_j\\ calibrated so the
effective degrees of freedom of the base-learner equals `df_target`.
Near-constant predictors fall back automatically to unpenalised linear
base-learners (with a warning).

`bgam()` complements the rest of the roadrunner lineup: it is more
interpretable than
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
(explicit additive components, no kernel matrix), handles non-linearity
without manual basis specification unlike
[`ols()`](https://cetialphafive.github.io/roadrunner/reference/ols.md),
and targets smooth monotonic / additive signals where
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)'s
piecewise-linear hinge functions may be less efficient. The complexity
per fit is O(n \* p \* K) per iteration (K ~ `nknots + degree + 1`),
making it practical for large-n, moderately large-p settings where
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
becomes memory-bound.

The formula method builds a standard model frame (factors expanded via
`model.matrix`, formula's `-1` disables the intercept initialisation),
stores `$terms` and `$xlevels` for
[`predict()`](https://rdrr.io/r/stats/predict.html), and dispatches to
the default method. The default method takes a numeric matrix directly.

## References

Buehlmann, P. and Yu, B. (2003). Boosting with the L2 loss: Regression
and classification. *Journal of the American Statistical Association*,
98(462):324–339.
[doi:10.1198/016214503000125](https://doi.org/10.1198/016214503000125)

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2):89–121.
[doi:10.1214/ss/1038425655](https://doi.org/10.1214/ss/1038425655)

Hofner, B., Mayr, A., Robinzonov, N. and Schmid, M. (2014). Model-based
boosting in R: A hands-on tutorial using the R package mboost.
*Computational Statistics*, 29(1–2):3–35.
[doi:10.1007/s00180-012-0382-5](https://doi.org/10.1007/s00180-012-0382-5)

Schmid, M. and Hothorn, T. (2008). Boosting additive models using
component-wise P-splines. *Computational Statistics and Data Analysis*,
53(2):298–311.
[doi:10.1016/j.csda.2008.09.009](https://doi.org/10.1016/j.csda.2008.09.009)

## Examples

``` r
set.seed(1)
n <- 200
x <- matrix(rnorm(n * 4), n, 4)
y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)

# Fit with fixed mstop (fast, no CV)
fit <- bgam(x, y, mstop = 100, autotune = FALSE, nthreads = 1)
print(fit)
#> Call:
#> bgam.default(x = x, y = y, mstop = 100, autotune = FALSE, nthreads = 1)
#> 
#> bgam fit:  family = gaussian 
#>   n = 200   p (predictors) = 4 
#>   mstop = 100   nu = 0.1 
#>   nknots = 20   degree = 3   dpen = 2   df_target = 4 
#> 
#> Top-5 predictors by selection frequency:
#>    V2           : 0.58 
#>    V1           : 0.27 
#>    V3           : 0.1 
#>    V4           : 0.05 
#> 
#> Residual sigma: 0.6854 
predict(fit, x[1:3, ])
#> [1] -0.3831642  1.6732858  0.7409848

# Formula interface
df <- data.frame(y = y, x)
fit2 <- bgam(y ~ ., data = df, mstop = 50, autotune = FALSE)
head(fit2$selection_frequency)
#>  X1  X2  X3  X4 
#> 0.4 0.6 0.0 0.0 
```
