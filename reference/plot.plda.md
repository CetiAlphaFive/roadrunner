# Projection plot for `plda` fits

Projection plot for `plda` fits

## Usage

``` r
# S3 method for class 'plda'
plot(x, data = NULL, labels = NULL, ...)
```

## Arguments

- x:

  A `"plda"` object.

- data:

  Predictor matrix (required).

- labels:

  Optional class label vector for colouring points.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
  plot(fit, data = iris[, 1:4])
} # }
```
