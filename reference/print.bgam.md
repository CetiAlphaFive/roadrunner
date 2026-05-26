# Print a `bgam` fit

Prints a compact summary of a fitted component-wise P-spline boosting
model: the call, family, sample size, number of base-learners, boosting
parameters (`mstop`, `nu`, `nknots`, `degree`, `dpen`, `df_target`),
CV-selected `mstop_opt` (when `autotune = TRUE`), bagging info (when
`n.boot > 0`), the top-5 predictors ranked by selection frequency, and
the residual sigma for gaussian fits.

## Usage

``` r
# S3 method for class 'bgam'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"bgam"`, as returned by
  [`bgam()`](https://cetialphafive.github.io/roadrunner/reference/bgam.md).

- digits:

  Number of significant digits for numeric output. Default
  `max(3L, getOption("digits") - 3L)`.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.

## See also

[`bgam()`](https://cetialphafive.github.io/roadrunner/reference/bgam.md),
[`summary.bgam()`](https://cetialphafive.github.io/roadrunner/reference/summary.bgam.md),
[`plot.bgam()`](https://cetialphafive.github.io/roadrunner/reference/plot.bgam.md)

## Examples

``` r
set.seed(1)
n <- 200
x <- matrix(rnorm(n * 4), n, 4)
y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
fit <- bgam(x, y, mstop = 50, autotune = FALSE)
print(fit)
#> Call:
#> bgam.default(x = x, y = y, mstop = 50, autotune = FALSE)
#> 
#> bgam fit:  family = gaussian 
#>   n = 200   p (predictors) = 4 
#>   mstop = 50   nu = 0.1 
#>   nknots = 20   degree = 3   dpen = 2   df_target = 4 
#> 
#> Top-5 predictors by selection frequency:
#>    V2           : 0.6 
#>    V1           : 0.4 
#>    V3           : 0 
#>    V4           : 0 
#> 
#> Residual sigma: 0.7089 
```
