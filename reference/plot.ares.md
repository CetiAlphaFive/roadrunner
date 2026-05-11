# Plot method for `ares` fits (residuals vs fitted)

Plot method for `ares` fits (residuals vs fitted)

## Usage

``` r
# S3 method for class 'ares'
plot(x, which = 1L, ...)
```

## Arguments

- x:

  an `ares` object

- which:

  integer in 1:1 — currently only residuals-vs-fitted is supported

- ...:

  passed to [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
  plot(fit)
} # }
```
