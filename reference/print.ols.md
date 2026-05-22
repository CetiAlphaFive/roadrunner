# Print method for `ols` fits

Prints the call, the fitted coefficients, the residual standard error,
and the residual degrees of freedom.

## Usage

``` r
# S3 method for class 'ols'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"ols"`.

- digits:

  Significant digits for numeric output.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
fit <- ols(mpg ~ wt + hp, data = mtcars)
print(fit)
#> Call:
#> ols.formula(x = mpg ~ wt + hp, data = mtcars)
#> 
#> linear model (OLS)
#> 
#> Coefficients:
#> (Intercept)           wt           hp  
#>    37.22727     -3.87783     -0.03177  
#> 
#>   Residual SE: 2.593 on 29 degrees of freedom
```
