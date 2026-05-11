# Summary method for `ares` fits

Summary method for `ares` fits

## Usage

``` r
# S3 method for class 'ares'
summary(object, ...)

# S3 method for class 'summary.ares'
print(x, ...)
```

## Arguments

- object:

  an `ares` object

- ...:

  unused

- x:

  a `summary.ares` object

## Value

An object of class `"summary.ares"`.

## Examples

``` r
fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
summary(fit)
#> Call:
#> ares.default(x = as.matrix(mtcars[, -1]), y = mtcars$mpg, nthreads = 2)
#> 
#> MARS summary:  family = gaussian 
#>   Selected: 6 of 21 forward-pass terms
#>   RSS = 91.81364   GCV = 6.662215   degree = 1 
#> 
#>           term        coef
#>    (Intercept) 19.71062781
#>    h(145-disp)  0.08426401
#>      h(123-hp)  0.08435673
#>      h(gear-4)  3.85511175
#>  h(17.05-qsec) -2.45008107
#>     h(wt-2.78) -2.95313879
```
