# Predict from a penalized LDA fit

Predict from a penalized LDA fit

## Usage

``` r
# S3 method for class 'plda'
predict(object, newdata, type = c("class", "posterior", "projection"), ...)
```

## Arguments

- object:

  A `"plda"` object.

- newdata:

  Numeric matrix or data frame of predictors. For formula-trained fits,
  `newdata` must be a data frame or a named matrix with the original
  predictor columns (matching the behavior of `predict.krls_rr`).

- type:

  `"class"` (default), `"posterior"`, or `"projection"`.

- ...:

  Unused.

## Value

Factor of class labels, posterior probability matrix, or projection
matrix (with columns named `"Discriminant 1"`, `"Discriminant 2"`, etc.)
when `type = "projection"`.

## Examples

``` r
fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
predict(fit, iris)
#>   [1] setosa     setosa     setosa     setosa     setosa     setosa    
#>   [7] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [13] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [19] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [25] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [31] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [37] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [43] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [49] setosa     setosa     versicolor versicolor versicolor versicolor
#>  [55] versicolor versicolor versicolor versicolor versicolor versicolor
#>  [61] versicolor versicolor versicolor versicolor versicolor versicolor
#>  [67] versicolor versicolor versicolor versicolor virginica  versicolor
#>  [73] versicolor versicolor versicolor versicolor versicolor virginica 
#>  [79] versicolor versicolor versicolor versicolor versicolor versicolor
#>  [85] versicolor versicolor versicolor versicolor versicolor versicolor
#>  [91] versicolor versicolor versicolor versicolor versicolor versicolor
#>  [97] versicolor versicolor versicolor versicolor virginica  virginica 
#> [103] virginica  virginica  virginica  virginica  versicolor virginica 
#> [109] virginica  virginica  virginica  virginica  virginica  virginica 
#> [115] virginica  virginica  virginica  virginica  virginica  versicolor
#> [121] virginica  virginica  virginica  virginica  virginica  virginica 
#> [127] virginica  virginica  virginica  virginica  virginica  virginica 
#> [133] virginica  versicolor versicolor virginica  virginica  virginica 
#> [139] virginica  virginica  virginica  virginica  virginica  virginica 
#> [145] virginica  virginica  virginica  virginica  virginica  virginica 
#> Levels: setosa versicolor virginica
```
