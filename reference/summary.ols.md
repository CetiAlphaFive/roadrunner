# Summary method for `ols` fits

Returns the coefficient table (estimate, standard error, t value, p
value), the residual standard error, R-squared and adjusted R-squared,
and the overall F-statistic.

## Usage

``` r
# S3 method for class 'ols'
summary(object, robust = c("none", "HC0", "HC1", "HC2", "HC3"), ...)

# S3 method for class 'summary.ols'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"ols"`.

- robust:

  Heteroscedasticity-consistent covariance for the standard errors and t
  / F tests. `"none"` (default) uses the classical covariance;
  `"HC0"`-`"HC3"` use a sandwich estimator.

- ...:

  Currently ignored.

- x:

  A `summary.ols` object.

- digits:

  Significant digits for numeric output.

## Value

An object of class `"summary.ols"`.

## Examples

``` r
fit <- ols(mpg ~ wt + hp, data = mtcars)
summary(fit)
#> Call:
#> ols.formula(x = mpg ~ wt + hp, data = mtcars)
#> 
#> linear model (OLS)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 37.22727    1.59879  23.285  < 2e-16 ***
#> wt          -3.87783    0.63273  -6.129 1.12e-06 ***
#> hp          -0.03177    0.00903  -3.519  0.00145 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 2.593 on 29 degrees of freedom
#> Multiple R-squared: 0.8268,  Adjusted R-squared: 0.8148
#> F-statistic: 69.21 on 2 and 29 DF,  p-value: 9.109e-12
summary(fit, robust = "HC3")
#> Call:
#> ols.formula(x = mpg ~ wt + hp, data = mtcars)
#> 
#> linear model (OLS)  --  HC3 robust standard errors
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 37.227270   2.229805  16.695  < 2e-16 ***
#> wt          -3.877831   0.768519  -5.046 2.23e-05 ***
#> hp          -0.031773   0.009385  -3.385  0.00206 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 2.593 on 29 degrees of freedom
#> Multiple R-squared: 0.8268,  Adjusted R-squared: 0.8148
#> F-statistic: 35.73 on 2 and 29 DF,  p-value: 1.499e-08
```
