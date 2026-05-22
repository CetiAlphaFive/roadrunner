# Print method for `logreg` fits

Prints the call, the fitted coefficients, the residual and null
deviance, the IRLS iteration count, and the convergence status.

## Usage

``` r
# S3 method for class 'logreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"logreg"`.

- digits:

  Significant digits for numeric output.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
print(logreg(y ~ wt + hp, data = df))
#> Call:
#> logreg.formula(x = y ~ wt + hp, data = df)
#> 
#> binary logistic regression
#> 
#> Coefficients:
#> (Intercept)           wt           hp  
#>    18.86630     -8.08348      0.03626  
#> 
#>   Null deviance:    43.23 on 31 degrees of freedom
#>   Residual deviance:10.06 on 29 degrees of freedom
#>   AIC: 16.06 
#>   IRLS iterations: 8 (converged) 
```
