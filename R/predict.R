# predict.ares -- bagging-aware prediction (gaussian + binomial)
#
# Three code paths:
#   1. newdata = NULL  -> return $fitted.values (training-set predictions).
#      With se.fit = TRUE and a bag, attaches a placeholder "sd" attribute
#      (bag SD on training rows requires the original x; see .ares_boot_sd).
#   2. newdata supplied, no $boot slot -> single-fit prediction via
#      mars_basis_cpp + coefficient multiply.
#   3. newdata supplied, $boot present -> average prediction across the
#      central fit plus each bootstrap replicate. With se.fit = TRUE, the
#      per-row sample SD across replicates is attached as attr(., "sd").
#
# When the fit was built with `family = "binomial"`, predictions go through
# the logit link. `type = "link"` returns the linear predictor x %*% beta;
# `type = "response"` returns plogis(linear predictor). For bag predictions
# under binomial, replicate-level predictions are averaged on the response
# (probability) scale so the bag mean is itself a valid probability.

#' Predictions from an `ares` fit
#'
#' @param object An object of class `"ares"`.
#' @param newdata A numeric matrix or data frame of new predictors. If `NULL`
#'   (default), returns `object$fitted.values`.
#' @param type Prediction scale. `"response"` (default) returns
#'   `plogis(eta)` for binomial fits and the fitted mean for gaussian fits;
#'   `"link"` returns the linear predictor `eta`. Ignored for gaussian
#'   models (link and response coincide).
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
predict.ares <- function(object, newdata = NULL,
                         type = c("response", "link"),
                         se.fit = FALSE, ...) {
  type <- match.arg(type)
  fam  <- if (is.null(object$family)) "gaussian" else object$family
  if (is.null(newdata)) {
    yhat <- if (fam == "binomial" && type == "link") {
      object$linear.predictor %||% object$fitted.values
    } else {
      object$fitted.values
    }
    if (isTRUE(se.fit) && !is.null(object$boot)) {
      # Re-derive bag predictions on the training data: each replicate's
      # fit holds dirs/cuts/selected.terms; we evaluate them via the basis
      # builder on the full design (the column space is set by namesx).
      attr(yhat, "sd") <- .ares_boot_sd(object, NULL)
    }
    return(yhat)
  }
  if (is.data.frame(newdata)) {
    # When the training fit expanded factor / character columns through
    # `model.matrix`, replay that expansion here using the stored
    # `$factor_info$xlevels`. The result is a numeric matrix whose
    # columns line up with the training design.
    if (!is.null(object$factor_info)) {
      fi <- object$factor_info
      # Require the original (pre-expansion) column names to be present.
      if (!all(fi$orig_names %in% colnames(newdata)))
        stop("ares: newdata is missing columns: ",
             paste(setdiff(fi$orig_names, colnames(newdata)), collapse = ", "))
      newdata <- newdata[, fi$orig_names, drop = FALSE]
      # Replay character/factor handling: coerce character to factor with
      # the *training* levels; refactor existing factor columns onto the
      # training levels. Out-of-vocabulary levels become NA -- they will
      # be picked up by the predict-side NA handling below.
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if (is.character(col) || is.factor(col)) {
          newdata[[jname]] <- factor(col, levels = fi$xlevels[[jname]])
        }
      }
      xnew <- stats::model.matrix(~ ., data = newdata)
      if ("(Intercept)" %in% colnames(xnew))
        xnew <- xnew[, colnames(xnew) != "(Intercept)", drop = FALSE]
      # Defensive: line up to training expansion exactly.
      if (!identical(colnames(xnew), fi$expanded_names)) {
        missing_cols <- setdiff(fi$expanded_names, colnames(xnew))
        if (length(missing_cols))
          stop("ares: newdata expansion is missing columns: ",
               paste(missing_cols, collapse = ", "),
               " -- likely an out-of-vocabulary factor level.")
        xnew <- xnew[, fi$expanded_names, drop = FALSE]
      }
    } else {
      if (!all(object$namesx %in% colnames(newdata)))
        stop("ares: newdata is missing columns: ",
             paste(setdiff(object$namesx, colnames(newdata)), collapse = ", "))
      xnew <- as.matrix(newdata[, object$namesx, drop = FALSE])
    }
  } else if (is.matrix(newdata) || is.numeric(newdata)) {
    xnew <- as.matrix(newdata)
    if (ncol(xnew) != length(object$namesx))
      stop("ares: newdata has ", ncol(xnew),
           " columns; expected ", length(object$namesx), ".")
  } else {
    stop("ares: newdata must be a matrix or data frame.")
  }
  storage.mode(xnew) <- "double"

  # Reapply the fit's NA-imputation to newdata so predictions don't go NA
  # downstream. The medians are the same ones used at fit time (stored in
  # `$na.medians` by ares.default when `na.action = "impute"`). If the
  # fit used `na.action = "omit"`, no medians were stored; in that case
  # newdata NAs propagate to NA predictions and we just warn.
  if (any(is.na(xnew))) {
    nrow_aff <- sum(rowSums(is.na(xnew)) > 0L)
    if (!is.null(object$na.medians)) {
      for (j in seq_len(ncol(xnew))) {
        na_j <- is.na(xnew[, j])
        if (any(na_j)) xnew[na_j, j] <- object$na.medians[j]
      }
      warning("Missing values in newdata: median-imputed ", nrow_aff,
              " row(s) using the column medians stored from training.",
              call. = FALSE)
    } else {
      warning("Missing values in newdata but no medians stored",
              " (training used na.action='omit'); the affected ",
              nrow_aff, " row(s) will return NA predictions.",
              call. = FALSE)
    }
  }

  bx_new <- mars_basis_cpp(xnew, object$dirs, object$cuts, object$selected.terms)
  eta_central <- drop(bx_new %*% object$coefficients)

  # Single-fit path.
  if (is.null(object$boot)) {
    if (fam == "binomial") {
      return(if (type == "link") eta_central else stats::plogis(eta_central))
    }
    return(eta_central)
  }

  # Bag predictions: average across (n.boot + 1) fits — the central fit plus
  # each replicate. For binomial we average on the response (probability)
  # scale so the bag mean stays in [0, 1]; for gaussian we average on the
  # raw scale (link == response).
  preds <- matrix(NA_real_, nrow = nrow(xnew),
                  ncol = length(object$boot$fits) + 1L)
  if (fam == "binomial") {
    preds[, 1] <- stats::plogis(eta_central)
    for (b in seq_along(object$boot$fits)) {
      fb <- object$boot$fits[[b]]
      bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
      eta_b <- drop(bx_b %*% fb$coefficients)
      preds[, b + 1L] <- stats::plogis(eta_b)
    }
  } else {
    preds[, 1] <- eta_central
    for (b in seq_along(object$boot$fits)) {
      fb <- object$boot$fits[[b]]
      bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
      preds[, b + 1L] <- drop(bx_b %*% fb$coefficients)
    }
  }
  yhat <- rowMeans(preds)
  if (isTRUE(se.fit)) {
    attr(yhat, "sd") <- apply(preds, 1, stats::sd)
  }
  # Convert to link scale only after averaging (for binomial): the
  # bag-mean on response scale is the natural prediction; if the user
  # explicitly asked type="link", invert via qlogis (handling 0/1 edge
  # cases by clamping to a small epsilon so the logit is finite).
  if (fam == "binomial" && type == "link") {
    eps <- 1e-15
    yhat <- stats::qlogis(pmin(pmax(yhat, eps), 1 - eps))
  }
  yhat
}

# Tiny null-coalescing operator (mirrors the one in ares.R).
`%||%` <- function(a, b) if (is.null(a)) b else a

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
