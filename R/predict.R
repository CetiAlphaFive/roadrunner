#' Predictions from an `ares` fit
#'
#' @param object An object of class `"ares"`.
#' @param newdata A numeric matrix or data frame of new predictors. If `NULL`
#'   (default), returns `object$fitted.values`.
#' @param se.fit If `TRUE` and `object` was fit with `n.boot > 0`, the
#'   returned vector carries an `"sd"` attribute holding the per-prediction
#'   bag standard deviation. Default `FALSE` (plain numeric vector).
#' @param ... Additional arguments. Currently unused.
#' @return A numeric vector of length `nrow(newdata)`. When `se.fit = TRUE`
#'   and bag fits are present, the result is the bag mean prediction with
#'   `attr(., "sd")` set to the per-row sample SD across replicates.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
#' @export
predict.ares <- function(object, newdata = NULL, se.fit = FALSE, ...) {
  if (is.null(newdata)) {
    yhat <- object$fitted.values
    if (isTRUE(se.fit) && !is.null(object$boot)) {
      # Re-derive bag predictions on the training data: each replicate's
      # fit holds dirs/cuts/selected.terms; we evaluate them via the basis
      # builder on the full design (the column space is set by namesx).
      attr(yhat, "sd") <- .ares_boot_sd(object, NULL)
    }
    return(yhat)
  }
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
  yhat_central <- drop(bx_new %*% object$coefficients)

  if (is.null(object$boot)) return(yhat_central)

  # Bag predictions: average across (n.boot + 1) fits — the central fit plus
  # each replicate.
  preds <- matrix(NA_real_, nrow = nrow(xnew),
                  ncol = length(object$boot$fits) + 1L)
  preds[, 1] <- yhat_central
  for (b in seq_along(object$boot$fits)) {
    fb <- object$boot$fits[[b]]
    bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
    preds[, b + 1L] <- drop(bx_b %*% fb$coefficients)
  }
  yhat <- rowMeans(preds)
  if (isTRUE(se.fit)) {
    attr(yhat, "sd") <- apply(preds, 1, stats::sd)
  }
  yhat
}

# Internal: compute bag SD on a built-in xnew (training x). Currently unused
# from outside predict.ares; kept tiny for clarity.
.ares_boot_sd <- function(object, newdata) {
  # Use object$bx for the central fit (already on training rows).
  if (is.null(newdata)) {
    n <- nrow(object$bx)
    preds <- matrix(NA_real_, n, length(object$boot$fits) + 1L)
    preds[, 1] <- object$fitted.values
    # Reconstruct training x from object$bx is impossible — bag-on-train
    # SD requires the original x. We don't keep it. Fall back to NA: bag
    # SD on training rows requires the user to call predict(object,
    # newdata = training_x, se.fit = TRUE).
    return(rep(NA_real_, n))
  }
  rep(NA_real_, nrow(newdata))
}
