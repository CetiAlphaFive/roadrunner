# Summary method for `logreg` fits

Returns the coefficient table (estimate, standard error, z value, p
value), the residual and null deviance, the AIC, and the IRLS iteration
count / convergence status.

## Usage

``` r
# S3 method for class 'logreg'
summary(object, robust = c("none", "HC0", "HC1", "HC2", "HC3"), ...)

# S3 method for class 'summary.logreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"logreg"`.

- robust:

  Heteroscedasticity-consistent covariance for the standard errors and z
  tests. `"none"` (default) uses the classical maximum-likelihood
  covariance; `"HC0"`-`"HC3"` use a sandwich estimator.

- ...:

  Currently ignored.

- x:

  A `summary.logreg` object.

- digits:

  Significant digits for numeric output.

## Value

An object of class `"summary.logreg"`.

## Examples

``` r
df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
fit <- logreg(y ~ wt + hp, data = df)
summary(fit)
#> Call:
#> logreg.formula(x = y ~ wt + hp, data = df)
#> 
#> binary logistic regression
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) 18.86630    7.44356   2.535  0.01126 * 
#> wt          -8.08348    3.06868  -2.634  0.00843 **
#> hp           0.03626    0.01773   2.044  0.04091 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>     Null deviance: 43.23  on 31 degrees of freedom
#> Residual deviance: 10.06  on 29 degrees of freedom
#> AIC: 16.06
#> Number of Fisher Scoring iterations: 8
summary(fit, robust = "HC3")
#> Call:
#> logreg.formula(x = y ~ wt + hp, data = df)
#> 
#> binary logistic regression  --  HC3 robust standard errors
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) 18.86630    9.76509   1.932  0.05336 . 
#> wt          -8.08348    3.29425  -2.454  0.01413 * 
#> hp           0.03626    0.01182   3.067  0.00216 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>     Null deviance: 43.23  on 31 degrees of freedom
#> Residual deviance: 10.06  on 29 degrees of freedom
#> AIC: 16.06
#> Number of Fisher Scoring iterations: 8
```
