# Summary method for `plda` fits

Summary method for `plda` fits

## Usage

``` r
# S3 method for class 'plda'
summary(object, ...)

# S3 method for class 'summary.plda'
print(x, ...)
```

## Arguments

- object:

  A `"plda"` object.

- ...:

  Currently ignored.

- x:

  A `summary.plda` object.

## Value

An object of class `"summary.plda"`.

Invisibly returns `x`.

## Examples

``` r
fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
summary(fit)
#> Penalized LDA summary
#>   penalty: L1 | classes: 3 | discriminants: 2 
#>   lambda: 0.1 
#>   discriminant 1: 4 / 4 nonzero features
#>   discriminant 2: 4 / 4 nonzero features
fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
print(summary(fit))
#> Penalized LDA summary
#>   penalty: L1 | classes: 3 | discriminants: 2 
#>   lambda: 0.1 
#>   discriminant 1: 4 / 4 nonzero features
#>   discriminant 2: 4 / 4 nonzero features
```
