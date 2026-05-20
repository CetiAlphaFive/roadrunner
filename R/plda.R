# R/plda.R

#' Penalized linear discriminant analysis
#'
#' Penalized Fisher's linear discriminant (Witten & Tibshirani 2011) with L1 or
#' fused-lasso penalties, multi-class support, and built-in cross-validation.
#'
#' @param x Numeric predictor matrix, or a formula for the formula interface.
#' @param y Factor (or coercible) class label of length `nrow(x)`.
#' @param ... Passed to methods.
#' @return An object of S3 class `"plda"`.
#' @export
plda <- function(x, ...) UseMethod("plda")
