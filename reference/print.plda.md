# Print method for `plda` fits

Print method for `plda` fits

## Usage

``` r
# S3 method for class 'plda'
print(x, ...)
```

## Arguments

- x:

  A `"plda"` object.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
print(fit)
#> Penalized LDA (L1 penalty)
#>   Classes:      setosa, versicolor, virginica
#>   Discriminants: 2 
#>   lambda:       0.1
#>   Nonzero feats per discriminant: 4, 4 
```
