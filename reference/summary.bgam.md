# Summarise a `bgam` fit

Returns a `summary.bgam` S3 object containing the call, family, sample
and predictor dimensions, boosting parameters (`mstop`, `nu`,
`df_target`), CV-selected `mstop_opt` (when `autotune = TRUE`), bagging
info (when `n.boot > 0`), training loss at `mstop`, residual sigma for
gaussian fits, and a full predictor selection-frequency table sorted in
descending order. The selection frequency is the fraction of boosting
iterations in which each predictor was the active base-learner.

## Usage

``` r
# S3 method for class 'bgam'
summary(object, ...)

# S3 method for class 'summary.bgam'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"bgam"`, as returned by
  [`bgam()`](https://cetialphafive.github.io/roadrunner/reference/bgam.md).

- ...:

  Currently ignored.

- x:

  A `summary.bgam` object returned by `summary.bgam()`.

- digits:

  Number of significant digits for numeric output. Default
  `max(3L, getOption("digits") - 3L)`.

## Value

An object of class `"summary.bgam"`. Print it with
`print.summary.bgam()`.

Invisibly returns `x`.

## See also

[`bgam()`](https://cetialphafive.github.io/roadrunner/reference/bgam.md),
[`print.bgam()`](https://cetialphafive.github.io/roadrunner/reference/print.bgam.md),
[`plot.bgam()`](https://cetialphafive.github.io/roadrunner/reference/plot.bgam.md)

## Examples

``` r
set.seed(1)
n <- 200
x <- matrix(rnorm(n * 4), n, 4)
y <- sin(x[, 1]) + 0.5 * x[, 2]^2 + rnorm(n, sd = 0.5)
fit <- bgam(x, y, mstop = 50, autotune = FALSE)
summary(fit)
#> Call:
#> bgam.default(x = x, y = y, mstop = 50, autotune = FALSE)
#> 
#> bgam summary:  family = gaussian 
#>   n = 200   p = 4 
#>   mstop = 50   nu = 0.1   df_target = 4 
#>   Training mse at mstop: 0.4925 
#>   Residual sigma: 0.7089 
#> 
#> Predictor selection frequencies (sorted):
#>  predictor sel_frequency
#>         V2           0.6
#>         V1           0.4
#>         V3           0.0
#>         V4           0.0
```
