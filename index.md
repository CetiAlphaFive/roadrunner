# roadrunner

Fast, low-dependency machine learning algorithms in R.

## Install

``` r

# install.packages("devtools")
devtools::install_github("CetiAlphaFive/roadrunner")
```

## Usage

``` r

library(roadrunner)

fit <- ares(mpg ~ ., data = mtcars, degree = 2)
predict(fit, head(mtcars))
```

## License

MIT (c) Jack Trametta
