# Diagnostic plots for an `ols` fit

Four diagnostic panels modelled on
[`stats::plot.lm()`](https://rdrr.io/r/stats/plot.lm.html): residuals vs
fitted, normal Q-Q of standardized residuals, scale-location, and
residuals vs leverage with Cook's-distance contours.

## Usage

``` r
# S3 method for class 'ols'
plot(x, which = 1:4, id.n = 3L, ...)
```

## Arguments

- x:

  A fitted object of class `"ols"`.

- which:

  Integer subset of `1:4` selecting which panels to draw. Default `1:4`.

- id.n:

  Number of extreme points to label per panel (default `3`; set `0` to
  suppress).

- ...:

  Further graphical parameters passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) calls.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  plot(fit)
} # }
```
