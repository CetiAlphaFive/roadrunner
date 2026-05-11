#' roadrunner: Fast, Low-Dependency Machine Learning Algorithms
#'
#' `roadrunner` is a collection of fast, low-dependency implementations of
#' classical machine learning algorithms. The first algorithm exposed is
#' [ares()], a Multivariate Adaptive Regression Splines (MARS) fitter built
#' on `Rcpp` and `RcppParallel`. Additional algorithms will be added as
#' separate exported functions under the same package roof.
#'
#' @keywords internal
#' @aliases roadrunner-package
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics*, 19(1):1-67.
"_PACKAGE"

## usethis namespace: start
#' @useDynLib roadrunner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel defaultNumThreads setThreadOptions
## usethis namespace: end
NULL
