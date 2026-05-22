# Diagnostic plots for a `logreg` fit

Four diagnostic panels for a logistic-regression fit: deviance residuals
vs the linear predictor, a normal Q-Q plot of the deviance residuals, a
scale-location panel, and deviance residuals vs leverage with
Cook's-distance contours.

## Usage

``` r
# S3 method for class 'logreg'
plot(x, which = 1:4, id.n = 3L, ...)
```

## Arguments

- x:

  A fitted object of class `"logreg"`.

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
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  plot(logreg(y ~ wt + hp, data = df))
} # }
```
