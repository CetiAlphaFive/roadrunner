#' Predictions from an `ares` fit
#'
#' Returns predictions for new data from a fitted `ares` MARS model.
#' For bagged fits (`n.boot > 0`), returns the mean across the central
#' fit and the bootstrap replicates; with `se.fit = TRUE`, the per-row
#' bag standard deviation is attached as an `"sd"` attribute. For
#' gaussian fits built with a residual variance model
#' (`varmod = "const"` or `"lm"`), `interval = "pint"` returns
#' prediction intervals.
#'
#' @param object An object of class `"ares"`.
#' @param newdata A numeric matrix or data frame of new predictors.
#'   `NULL` (default) returns `object$fitted.values`.
#' @param type Prediction scale.
#'   - `"response"` (default): probabilities for binomial, response-scale
#'     means for poisson and gamma, fitted values for gaussian.
#'   - `"link"`: the linear predictor. Coincides with `"response"` for
#'     gaussian fits.
#' @param se.fit If `TRUE` and the fit was built with `n.boot > 0`, the
#'   returned vector carries an `"sd"` attribute holding the per-row
#'   bag standard deviation. Default `FALSE`.
#' @param interval Type of interval to return.
#'   - `"none"` (default): return point predictions only.
#'   - `"pint"`: return a matrix with columns `c("fit", "lwr", "upr")`
#'     using the variance model stored at fit time. Requires
#'     `family = "gaussian"` and `varmod = "const"` or `"lm"`; errors
#'     otherwise.
#' @param level Confidence level for prediction intervals when
#'   `interval = "pint"`. Default `0.95`.
#' @param ... Currently unused.
#' @return A numeric vector of predictions, or a matrix with columns
#'   `c("fit", "lwr", "upr")` when `interval = "pint"`. With
#'   `se.fit = TRUE` on a bagged fit, the vector carries an `"sd"`
#'   attribute.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
#' @export
predict.ares <- function(object, newdata = NULL,
                         type = c("response", "link"),
                         se.fit = FALSE,
                         interval = c("none", "pint"),
                         level = 0.95, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  fam  <- if (is.null(object$family)) "gaussian" else object$family
  if (is.null(newdata)) {
    # BUG-002 (v0.0.0.9029): for bagged fits, `predict(object)` must return
    # the bag mean (matching `predict(object, x_train)`). Previously the
    # NULL branch short-circuited to `$fitted.values` (the central fit
    # only), silently diverging from the newdata-supplied path. We now
    # re-route through the bag path using the stored training x. The
    # short-circuit is kept for single (non-bagged) fits where central ==
    # bag-mean by definition.
    if (!is.null(object$boot) && !is.null(object$x)) {
      return(Recall(object, newdata = object$x, type = type,
                    se.fit = se.fit, interval = interval, level = level))
    }
    yhat <- if (fam == "binomial" && type == "link") {
      object$linear.predictor %||% object$fitted.values
    } else {
      object$fitted.values
    }
    if (isTRUE(se.fit) && !is.null(object$boot)) {
      attr(yhat, "sd") <- .ares_boot_sd(object, NULL)
    }
    if (interval == "pint") {
      if (fam != "gaussian")
        stop("ares: interval = 'pint' is only supported for",
             " family = 'gaussian' (got '", fam, "').")
      pi_mat <- .ares_make_pi(yhat, object, level)
      if (is.null(pi_mat))
        stop("ares: interval = 'pint' requires the fit was built",
             " with varmod = 'const' or 'lm'. Refit ares() with",
             " varmod = 'const' (or 'lm') and try again.")
      return(pi_mat)
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
      # Pre-expansion numeric-column NA fill using the training medians
      # captured in $factor_info$num_medians. Without this, the default
      # `model.matrix(~ ., newdata)` drops rows with any NA in numeric
      # columns and the result misaligns with the training expansion.
      if (!is.null(fi$num_medians)) {
        for (jn in names(fi$num_medians)) {
          if (jn %in% names(newdata)) {
            col <- newdata[[jn]]
            if (is.numeric(col) && anyNA(col)) {
              col[is.na(col)] <- fi$num_medians[[jn]]
              newdata[[jn]] <- col
            }
          }
        }
      }
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
    yhat <- if (fam == "binomial") {
      if (type == "link") eta_central
      else stats::plogis(pmin(pmax(eta_central, -30), 30))
    } else if (fam == "poisson" || fam == "gamma") {
      if (type == "link") eta_central
      else {
        clamp_hi <- if (fam == "poisson") 50 else 20
        clamp_lo <- if (fam == "poisson") -50 else -20
        exp(pmin(pmax(eta_central, clamp_lo), clamp_hi))
      }
    } else {
      eta_central
    }
    if (interval == "pint") {
      if (fam != "gaussian")
        stop("ares: interval = 'pint' is only supported for",
             " family = 'gaussian' (got '", fam, "').")
      pi_mat <- .ares_make_pi(yhat, object, level)
      if (is.null(pi_mat))
        stop("ares: interval = 'pint' requires the fit was built",
             " with varmod = 'const' or 'lm'.")
      return(pi_mat)
    }
    return(yhat)
  }

  # Bag predictions: average across (n.boot + 1) fits — the central fit plus
  # each replicate. For binomial we average on the response (probability)
  # scale so the bag mean stays in [0, 1]; for gaussian we average on the
  # raw scale (link == response).
  preds <- matrix(NA_real_, nrow = nrow(xnew),
                  ncol = length(object$boot$fits) + 1L)
  if (fam == "binomial") {
    preds[, 1] <- stats::plogis(pmin(pmax(eta_central, -30), 30))
    for (b in seq_along(object$boot$fits)) {
      fb <- object$boot$fits[[b]]
      bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
      eta_b <- drop(bx_b %*% fb$coefficients)
      preds[, b + 1L] <- stats::plogis(pmin(pmax(eta_b, -30), 30))
    }
  } else if (fam == "poisson" || fam == "gamma") {
    # Bag-average on the response (rate / mean) scale so the average stays
    # positive; convert to link scale (log of the mean) only if the user
    # asked for it. Matches the binomial probability-scale convention.
    clamp_hi <- if (fam == "poisson") 50 else 20
    clamp_lo <- if (fam == "poisson") -50 else -20
    preds[, 1] <- exp(pmin(pmax(eta_central, clamp_lo), clamp_hi))
    for (b in seq_along(object$boot$fits)) {
      fb <- object$boot$fits[[b]]
      bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
      eta_b <- drop(bx_b %*% fb$coefficients)
      preds[, b + 1L] <- exp(pmin(pmax(eta_b, clamp_lo), clamp_hi))
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
  # Convert to link scale only after averaging:
  #   - binomial: invert via qlogis(); clamp to (eps, 1-eps) for finiteness.
  #   - poisson/gamma: log(); clamp to a small positive epsilon so log is finite.
  if (fam == "binomial" && type == "link") {
    eps <- 1e-15
    yhat <- stats::qlogis(pmin(pmax(yhat, eps), 1 - eps))
  } else if ((fam == "poisson" || fam == "gamma") && type == "link") {
    eps <- 1e-300
    yhat <- log(pmax(yhat, eps))
  }
  if (interval == "pint") {
    if (fam != "gaussian")
      stop("ares: interval = 'pint' is only supported for",
           " family = 'gaussian' (got '", fam, "').")
    pi_mat <- .ares_make_pi(yhat, object, level)
    if (is.null(pi_mat))
      stop("ares: interval = 'pint' requires the fit was built",
           " with varmod = 'const' or 'lm'.")
    # PIs use the bag mean for `fit`; bag SD is approximate for the
    # variance model component but we don't combine them here (the
    # variance model already captures residual uncertainty). attr 'sd'
    # from bagging is dropped on the PI matrix path for clarity.
    return(pi_mat)
  }
  yhat
}

# Tiny null-coalescing operator (mirrors the one in ares.R).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Internal helper: turn a vector of fitted means into a (fit, lwr, upr)
# matrix using the variance model stored on `object`. Returns NULL when the
# variance model is missing, with the caller responsible for stop()ing.
.ares_make_pi <- function(yhat, object, level = 0.95) {
  vm <- object$varmod
  if (is.null(vm)) return(NULL)
  alpha <- (1 - level) / 2
  qq <- stats::qt(1 - alpha, df = vm$df)
  sigma_vec <- switch(vm$type,
    const = rep(vm$sigma_hat, length(yhat)),
    lm    = pmax(vm$scale * (vm$intercept + vm$slope * yhat),
                 1e-12)              # floor to keep sigma positive
  )
  cbind(fit = yhat, lwr = yhat - qq * sigma_vec,
        upr = yhat + qq * sigma_vec)
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
