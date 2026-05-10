#' ares: Fast Multivariate Adaptive Regression Splines
#'
#' `ares` fits MARS models using Friedman's (1991) fast least-squares forward
#' selection and a parallelized knot search built on `Rcpp` and `RcppParallel`.
#' Designed for numerical parity with the `earth` package on the gaussian-only
#' core while taking advantage of multi-core CPUs.
#'
#' @keywords internal
#' @aliases ares-package
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics*, 19(1):1-67.
"_PACKAGE"

## usethis namespace: start
#' @useDynLib ares, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel defaultNumThreads setThreadOptions
## usethis namespace: end
NULL
