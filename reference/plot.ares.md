# Plot method for `ares` fits

Diagnostic plot for a fitted `ares` model. Currently produces a
residuals-vs-fitted scatter; additional diagnostics may be added in the
future.

## Usage

``` r
# S3 method for class 'ares'
plot(x, which = 1L, ...)
```

## Arguments

- x:

  An object of class `"ares"`.

- which:

  Which diagnostic to draw. Only `1` (residuals vs fitted) is
  implemented.

- ...:

  Further graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
  plot(fit)
} # }
```
