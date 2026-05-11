# Print method for `ares` fits

Prints the call, family, term count, training RSS and GCV, and the
control parameters used. For `family = "binomial"`, also prints the
deviance, McFadden pseudo-R^2, AIC, and IRLS convergence.

## Usage

``` r
# S3 method for class 'ares'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"ares"`.

- digits:

  Significant digits for numeric output.

- ...:

  Currently ignored.

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
