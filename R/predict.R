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
      # BUG-010 (v0.0.0.9032): detect NA in factor / character columns up
      # front. Otherwise model.matrix(~ ., newdata) with na.action=na.omit
      # silently drops those rows and the result has fewer rows than
      # nrow(newdata) -- the predict()-length contract gets violated. We
      # error with a clear message naming the column(s), mirroring the
      # BUG-004 OOV pattern. (Numeric NAs are handled upstream via
      # num_medians or via the na.action="omit" mask in the matrix
      # branch below.)
      fna_report <- character(0)
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
                          sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                                  jname, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("ares: NA value(s) in factor/character newdata column(s); ",
             "cannot expand to the training design matrix. Drop the ",
             "offending row(s) or impute the categorical level(s) before ",
             "predicting.\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      # Replay character/factor handling: coerce character to factor with
      # the *training* levels; refactor existing factor columns onto the
      # training levels. BUG-004 (v0.0.0.9029): out-of-vocabulary levels
      # used to become NA and then model.matrix(..., na.action=na.omit)
      # silently dropped those ROWS from the result, so predict()
      # returned fewer rows than nrow(newdata). We now detect OOV per
      # column up front, collect offending row indices and level names,
      # and stop() with a clear error -- safer than silent row drops.
      oov_report <- list()
      for (jname in names(fi$xlevels)) {
        col <- newdata[[jname]]
        if (is.character(col) || is.factor(col)) {
          col_chr <- as.character(col)
          tr_lev <- fi$xlevels[[jname]]
          bad_rows <- which(!is.na(col_chr) & !(col_chr %in% tr_lev))
          if (length(bad_rows)) {
            bad_lev <- unique(col_chr[bad_rows])
            oov_report[[jname]] <- list(rows = bad_rows, levels = bad_lev)
          }
          newdata[[jname]] <- factor(col, levels = tr_lev)
        }
      }
      if (length(oov_report)) {
        msgs <- vapply(names(oov_report), function(jn) {
          rep <- oov_report[[jn]]
          sprintf("  column %s: %d row(s) (e.g. row %d) with level(s) not seen at fit time: %s",
                  jn, length(rep$rows), rep$rows[1],
                  paste(utils::head(rep$levels, 5), collapse = ", "))
        }, character(1L))
        stop("ares: out-of-vocabulary factor level(s) in newdata; ",
             "cannot expand to the training design matrix. ",
             "Fix the offending column(s) (e.g. drop those rows or re-fit ",
             "with the new level(s)) and retry.\n",
             paste(msgs, collapse = "\n"), call. = FALSE)
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
    } else if (!is.null(object$terms)) {
      # BUG-012 (v0.0.0.9032): formula-path fit with derived terms (I(x^2),
      # poly(x, 2), log(x + 10), scale(x), splines::bs(x), ...) needs the
      # original terms object to re-evaluate the design at predict time.
      # The previous fallback (line up newdata columns to object$namesx)
      # only worked when namesx referenced plain data-frame columns; it
      # errored on any derived term whose name didn't match a column of
      # newdata. We now drop the response side of the terms and rebuild
      # the design via model.matrix(..., xlev = object$xlevels) so factor
      # and character handling matches the training fit.
      tt <- stats::delete.response(object$terms)
      # Detect NA in factor columns referenced by the terms object before
      # model.matrix's default na.omit silently drops rows.
      term_vars <- all.vars(tt)
      fna_report <- character(0)
      for (jn in intersect(term_vars, colnames(newdata))) {
        col <- newdata[[jn]]
        if ((is.character(col) || is.factor(col)) && anyNA(col)) {
          bad_rows <- which(is.na(col))
          fna_report <- c(fna_report,
                          sprintf("  column %s: %d NA row(s) (e.g. row %d)",
                                  jn, length(bad_rows), bad_rows[1]))
        }
      }
      if (length(fna_report)) {
        stop("ares: NA value(s) in factor/character newdata column(s) ",
             "referenced by the model formula; cannot build the design ",
             "matrix. Drop the offending row(s) or impute the ",
             "categorical level(s) before predicting.\n",
             paste(fna_report, collapse = "\n"), call. = FALSE)
      }
      xlev <- if (!is.null(object$xlevels)) object$xlevels else NULL
      mm <- stats::model.matrix(tt, data = newdata, xlev = xlev)
      if ("(Intercept)" %in% colnames(mm))
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      # If the formula references columns we don't have, model.matrix
      # already errors; if it produced an unexpected column set, surface
      # a clean message and (defensively) re-order to the training set.
      if (!identical(colnames(mm), object$namesx)) {
        missing_cols <- setdiff(object$namesx, colnames(mm))
        if (length(missing_cols))
          stop("ares: newdata model-matrix expansion is missing columns: ",
               paste(missing_cols, collapse = ", "),
               " -- check that newdata has the same predictors used at fit time.")
        mm <- mm[, object$namesx, drop = FALSE]
      }
      xnew <- mm
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
  if (identical(vm$type, "const")) {
    sigma_vec <- rep(vm$sigma_hat, length(yhat))
  } else {
    # vm$type == "lm" -- linear MAD model in yhat. BUG-003 (v0.0.0.9029):
    # at extrapolation rows the predicted MAD can go non-positive, and
    # the prior code floored at 1e-12 -- which produced PI widths of
    # ~4e-12 with NO warning. We now floor at a meaningful lower bound
    # captured at fit time (max of 5% of the const-equivalent sigma and
    # half the smallest in-sample positive MAD), and warn when ANY row
    # hits that floor so the user knows the variance model is being
    # extrapolated past its validity range.
    raw_mad <- vm$scale * (vm$intercept + vm$slope * yhat)
    floor_val <- if (!is.null(vm$sigma_floor)) vm$sigma_floor else 1e-12
    needs_floor <- raw_mad < floor_val
    if (any(needs_floor)) {
      yhat_range_msg <- if (!is.null(vm$yhat_min) && !is.null(vm$yhat_max))
        sprintf(" (training yhat range [%.3g, %.3g])",
                vm$yhat_min, vm$yhat_max) else ""
      warning("ares: varmod = 'lm' produced non-positive predicted",
              " MAD at ", sum(needs_floor), " of ", length(yhat),
              " row(s)", yhat_range_msg,
              "; flooring PI sigma at ", signif(floor_val, 3),
              ". Predicted intervals at these rows reflect the floor",
              " (the variance model is being extrapolated outside its",
              " in-sample range).", call. = FALSE)
    }
    sigma_vec <- pmax(raw_mad, floor_val)
  }
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
