# Print method for `ares` fits

Print method for `ares` fits

## Usage

``` r
# S3 method for class 'ares'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an `ares` object

- digits:

  significant digits for numeric output

- ...:

  unused

## Value

Invisibly returns `x`.

## Examples

``` r
fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
print(fit)
#> Call:
#> ares.default(x = as.matrix(mtcars[, -1]), y = mtcars$mpg, nthreads = 2)
#> 
#> MARS fit:  family = gaussian 
#>   Selected terms: 6 of 21 forward-pass terms
#>   RSS: 91.81 
#>   GCV: 6.662 
#>   Degree: 1   Penalty: 2   nthreads: 2 
```
