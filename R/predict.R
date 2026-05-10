#' Predictions from an `ares` fit
#'
#' @param object An object of class `"ares"`.
#' @param newdata A numeric matrix or data frame of new predictors. If `NULL`
#'   (default), returns `object$fitted.values`.
#' @param ... Additional arguments. Currently unused.
#' @return A numeric vector of length `nrow(newdata)`.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
#' @export
predict.ares <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  if (is.data.frame(newdata)) {
    if (!all(object$namesx %in% colnames(newdata)))
      stop("ares: newdata is missing columns: ",
           paste(setdiff(object$namesx, colnames(newdata)), collapse = ", "))
    xnew <- as.matrix(newdata[, object$namesx, drop = FALSE])
  } else if (is.matrix(newdata) || is.numeric(newdata)) {
    xnew <- as.matrix(newdata)
    if (ncol(xnew) != length(object$namesx))
      stop("ares: newdata has ", ncol(xnew),
           " columns; expected ", length(object$namesx), ".")
  } else {
    stop("ares: newdata must be a matrix or data frame.")
  }
  storage.mode(xnew) <- "double"
  bx_new <- mars_basis_cpp(xnew, object$dirs, object$cuts, object$selected.terms)
  drop(bx_new %*% object$coefficients)
}
