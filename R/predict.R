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
      # BUG-015 (v0.0.0.9055): build the model frame with
      # na.action = na.pass so NAs in numeric columns flow through to
      # xnew and are handled downstream (median-imputation or post-basis
      # NA mask). model.matrix(<formula>, data = <df>) goes through
      # model.frame with the global na.action option (default na.omit),
      # which silently drops rows with NA -- breaking the
      # length(predict()) == nrow(newdata) contract.
      mf_new <- stats::model.frame(~ ., data = newdata,
                                   na.action = stats::na.pass)
      xnew <- stats::model.matrix(~ ., data = mf_new)
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
      # BUG-015 (v0.0.0.9055): build the model frame with na.pass so NAs
      # in numeric columns are preserved (the default na.omit silently
      # drops rows, breaking the predict()-length contract).
      mf_t <- stats::model.frame(tt, data = newdata, xlev = xlev,
                                 na.action = stats::na.pass)
      mm <- stats::model.matrix(tt, data = mf_t, xlev = xlev)
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
  # we honour the prior warning's promise and return NA predictions for
  # the affected rows (BUG-008, v0.0.0.9032 -- previously the warning
  # promised NA predictions but `NaN > 0` in the C++ hinge evaluates to
  # FALSE, silently collapsing the hinge to 0 and yielding finite WRONG
  # predictions for those rows). We track the affected row mask BEFORE
  # the C++ basis call, fill the offending cells with 0 so the C++ pass
  # doesn't propagate NaN/Inf into the basis, and re-impose NA on those
  # rows immediately after the linear-predictor compute below.
  na_rows_mask <- logical(nrow(xnew))
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
      # na.action = "omit" path: no medians stored. Remember the rows so
      # we can NA them out post-basis-evaluation, and zero-fill the
      # source cells so the C++ engine sees finite input.
      na_rows_mask <- rowSums(is.na(xnew)) > 0L
      xnew[is.na(xnew)] <- 0
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
    # BUG-008 (v0.0.0.9032): honour the omit-path NA promise.
    if (any(na_rows_mask)) yhat[na_rows_mask] <- NA_real_
    if (interval == "pint") {
      if (fam != "gaussian")
        stop("ares: interval = 'pint' is only supported for",
             " family = 'gaussian' (got '", fam, "').")
      pi_mat <- .ares_make_pi(yhat, object, level)
      if (is.null(pi_mat))
        stop("ares: interval = 'pint' requires the fit was built",
             " with varmod = 'const' or 'lm'.")
      if (any(na_rows_mask)) pi_mat[na_rows_mask, ] <- NA_real_
      return(pi_mat)
    }
    return(yhat)
  }

  # Bag predictions: collect per-replicate linear predictors `etas` AND
  # response-scale predictions `resps` separately. BUG-009 (v0.0.0.9032):
  # the link-scale return must be mean(eta_b), NOT g(mean(g^{-1}(eta_b))).
  # The latter (Jensen-violating) form used to drift by ~hundreds of
  # log-odds units on moderate-signal binomial bags, because
  # qlogis(mean(plogis(eta_b))) != mean(eta_b) whenever eta_b varies.
  # Response-scale return stays as before -- mean of the inverse-link
  # mapping is the right point estimate for type="response".
  nb1 <- length(object$boot$fits) + 1L
  etas  <- matrix(NA_real_, nrow = nrow(xnew), ncol = nb1)
  resps <- matrix(NA_real_, nrow = nrow(xnew), ncol = nb1)
  etas[, 1]  <- eta_central
  resps[, 1] <- if (fam == "binomial") {
    stats::plogis(pmin(pmax(eta_central, -30), 30))
  } else if (fam == "poisson" || fam == "gamma") {
    clamp_hi <- if (fam == "poisson") 50 else 20
    clamp_lo <- if (fam == "poisson") -50 else -20
    exp(pmin(pmax(eta_central, clamp_lo), clamp_hi))
  } else {
    eta_central
  }
  for (b in seq_along(object$boot$fits)) {
    fb <- object$boot$fits[[b]]
    bx_b <- mars_basis_cpp(xnew, fb$dirs, fb$cuts, fb$selected.terms)
    eta_b <- drop(bx_b %*% fb$coefficients)
    etas[, b + 1L] <- eta_b
    resps[, b + 1L] <- if (fam == "binomial") {
      stats::plogis(pmin(pmax(eta_b, -30), 30))
    } else if (fam == "poisson" || fam == "gamma") {
      clamp_hi <- if (fam == "poisson") 50 else 20
      clamp_lo <- if (fam == "poisson") -50 else -20
      exp(pmin(pmax(eta_b, clamp_lo), clamp_hi))
    } else {
      eta_b
    }
  }
  # Aggregate per `type`. Link scale -> mean of linear predictors (correct
  # Jensen-safe point estimate). Response scale -> mean of inverse-link
  # mappings (matches single-fit response). Gaussian: link == response.
  if (type == "link") {
    yhat <- rowMeans(etas)
    if (isTRUE(se.fit)) attr(yhat, "sd") <- apply(etas, 1, stats::sd)
  } else {
    yhat <- rowMeans(resps)
    if (isTRUE(se.fit)) attr(yhat, "sd") <- apply(resps, 1, stats::sd)
  }
  # BUG-008 (v0.0.0.9032): honour the omit-path NA promise on bag mean + SE.
  if (any(na_rows_mask)) {
    yhat[na_rows_mask] <- NA_real_
    if (isTRUE(se.fit)) {
      sdv <- attr(yhat, "sd")
      sdv[na_rows_mask] <- NA_real_
      attr(yhat, "sd") <- sdv
    }
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
    if (any(na_rows_mask)) pi_mat[na_rows_mask, ] <- NA_real_
    return(pi_mat)
  }
  yhat
}

#' Predict from a penalized LDA fit
#'
#' @param object A `"plda"` object.
#' @param newdata Numeric matrix or data frame of predictors. For
#'   formula-trained fits, `newdata` must be a data frame or a named
#'   matrix with the original predictor columns (matching the behavior
#'   of `predict.krls_rr`).
#' @param type `"class"` (default), `"posterior"`, or `"projection"`.
#' @param ... Unused.
#' @return Factor of class labels, posterior probability matrix, or
#'   projection matrix (with columns named `"Discriminant 1"`,
#'   `"Discriminant 2"`, etc.) when `type = "projection"`.
#' @examples
#' fit <- plda(Species ~ ., data = iris, lambda = 0.1, autotune = FALSE)
#' predict(fit, iris)
#' @export
predict.plda <- function(object, newdata,
                         type = c("class", "posterior", "projection"), ...) {
  type <- match.arg(type)
  if (!is.null(object$terms)) {
    trms <- stats::delete.response(object$terms)
    mf <- stats::model.frame(trms, as.data.frame(newdata),
                             xlev = object$xlevels)
    x <- stats::model.matrix(trms, mf)
    icpt <- match("(Intercept)", colnames(x), nomatch = 0L)
    if (icpt > 0L) x <- x[, -icpt, drop = FALSE]
  } else {
    x <- as.matrix(newdata)
    if (ncol(x) != nrow(object$discrim))
      stop(sprintf("plda: newdata has %d column(s); model expects %d.",
                   ncol(x), nrow(object$discrim)), call. = FALSE)
  }
  scores  <- plda_project_cpp(x, object$mu, object$sdw, object$discrim)
  colnames(scores) <- paste0("Discriminant ", seq_len(ncol(scores)))
  cscores <- object$cmeans %*% object$discrim   # G x K centroid projections
  if (type == "projection") return(scores)
  d2 <- outer(rowSums(scores^2), rowSums(cscores^2), `+`) -
    2 * scores %*% t(cscores)
  if (type == "posterior") {
    m <- d2 - apply(d2, 1L, min)
    e <- exp(-0.5 * m)
    p <- e / rowSums(e)
    colnames(p) <- object$classes
    return(p)
  }
  factor(object$classes[max.col(-d2, ties.method = "first")],
         levels = object$classes)
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

#' Predictions from an `ols` fit
#'
#' Returns predictions for new data from a fitted `ols` model, with
#' optional standard errors and exact closed-form confidence or
#' prediction intervals.
#'
#' @param object An object of class `"ols"`.
#' @param newdata A numeric matrix or data frame of new predictors.
#'   `NULL` (default) returns `object$fitted.values`.
#' @param type Prediction scale. Only `"response"` is supported (the
#'   linear predictor coincides with the response for a gaussian fit).
#' @param se.fit If `TRUE`, the returned vector carries an `"se.fit"`
#'   attribute holding the per-row standard error of the fitted mean.
#'   For bagged fits the standard error is the per-row bag standard
#'   deviation. Default `FALSE`.
#' @param interval Type of interval to return.
#'   - `"none"` (default): point predictions only.
#'   - `"confidence"`: a matrix with columns `c("fit", "lwr", "upr")`
#'     for the mean response.
#'   - `"prediction"`: a matrix with columns `c("fit", "lwr", "upr")`
#'     for a new observation (adds the residual variance, or the
#'     variance-model variance when the fit was built with
#'     `varmod = "const"` or `"lm"`).
#' @param level Confidence/prediction level. Default `0.95`.
#' @param robust Heteroscedasticity-consistent covariance used for the
#'   standard errors / intervals.
#'   - `"none"` (default): the classical `sigma2 (X'WX)^-1`.
#'   - `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`: a sandwich estimator.
#' @param ... Currently unused.
#' @return A numeric vector of predictions, or a matrix with columns
#'   `c("fit", "lwr", "upr")` when `interval != "none"`. With
#'   `se.fit = TRUE` the result carries an `"se.fit"` attribute.
#' @examples
#' fit <- ols(mpg ~ wt + hp, data = mtcars)
#' predict(fit, mtcars[1:5, ], interval = "confidence")
#' @export
predict.ols <- function(object, newdata = NULL,
                        type = "response",
                        se.fit = FALSE,
                        interval = c("none", "confidence", "prediction"),
                        level = 0.95,
                        robust = c("none", "HC0", "HC1", "HC2", "HC3"),
                        ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  robust <- match.arg(robust)

  # ---- newdata = NULL fast path -------------------------------------------
  # For a single (non-bagged) fit, predict(object) is just $fitted.values.
  # For a bagged fit it must be the bag mean (matching predict(object,
  # x_train)); we re-route through the bag path using the stored design.
  if (is.null(newdata)) {
    if (is.null(object$boot)) {
      yhat <- object$fitted.values
      if (!isTRUE(se.fit) && interval == "none") return(yhat)
      xnew <- object$X
    } else {
      xnew <- object$X
    }
  } else {
    xnew <- .ols_newdata_matrix(object, newdata)
  }
  storage.mode(xnew) <- "double"

  # ---- Point estimate -----------------------------------------------------
  # Single fit: the central-fit linear predictor. Bagged fit: the mean over
  # the central fit plus the bootstrap replicates (the bag mean).
  eta_central <- drop(xnew %*% object$coefficients)
  boot_preds <- NULL
  if (!is.null(object$boot)) {
    boot_preds <- cbind(eta_central, xnew %*% object$boot$coefficients)
    yhat <- rowMeans(boot_preds)
  } else {
    yhat <- eta_central
  }

  # ---- Variance-covariance of the coefficients ----------------------------
  # robust = "none" reuses the stored classical vcov; HC* recomputes the
  # sandwich on the training design via the C++ engine.
  vcov <- .ols_vcov(object, robust)

  need_var <- isTRUE(se.fit) || interval != "none"
  se_mean <- NULL
  if (need_var) {
    if (!is.null(object$boot)) {
      # Bagged fit: the per-row SE is the bag standard deviation of the
      # replicate predictions (plus the central fit).
      se_mean <- apply(boot_preds, 1L, stats::sd)
    } else {
      # var(xhat' beta) = xhat' V xhat, computed rowwise.
      se_mean <- sqrt(pmax(rowSums((xnew %*% vcov) * xnew), 0))
    }
  }

  if (isTRUE(se.fit)) attr(yhat, "se.fit") <- se_mean
  if (interval == "none") return(yhat)

  # ---- Closed-form confidence / prediction intervals ----------------------
  alpha <- (1 - level) / 2
  df <- object$df.residual
  qq <- stats::qt(1 - alpha, df = df)
  if (interval == "confidence") {
    half <- qq * se_mean
  } else {
    # Prediction interval: add the variance of a fresh observation. With a
    # variance model, use its (possibly yhat-dependent) sigma; otherwise
    # use the fitted residual variance sigma2.
    if (!is.null(object$varmod)) {
      sig_obs <- .ols_varmod_sigma(object$varmod, yhat)
      half <- qq * sqrt(se_mean^2 + sig_obs^2)
    } else {
      half <- qq * sqrt(se_mean^2 + object$sigma2)
    }
  }
  out <- cbind(fit = yhat, lwr = yhat - half, upr = yhat + half)
  if (isTRUE(se.fit)) attr(out, "se.fit") <- se_mean
  out
}

# Internal: turn `newdata` into a numeric design matrix aligned with the
# training design of an `ols` fit. Handles the formula path (terms +
# xlevels), the factor-expanded data-frame path, and the plain matrix /
# data-frame path. Mirrors the predict.ares column-alignment logic.
#
# @keywords internal
.ols_newdata_matrix <- function(object, newdata) {
  cn <- colnames(object$X)
  has_int <- isTRUE(object$intercept) ||
    (!is.null(cn) && "(Intercept)" %in% cn)
  if (!is.null(object$terms)) {
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, as.data.frame(newdata),
                             xlev = object$xlevels)
    mm <- stats::model.matrix(tt, mf)
    if (!identical(colnames(mm), cn)) {
      missing_cols <- setdiff(cn, colnames(mm))
      if (length(missing_cols))
        stop("ols: newdata model-matrix expansion is missing columns: ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
      mm <- mm[, cn, drop = FALSE]
    }
    return(mm)
  }
  if (is.data.frame(newdata) && !is.null(object$factor_info)) {
    fi <- object$factor_info
    if (!all(fi$orig_names %in% colnames(newdata)))
      stop("ols: newdata is missing columns: ",
           paste(setdiff(fi$orig_names, colnames(newdata)),
                 collapse = ", "), call. = FALSE)
    nd <- newdata[, fi$orig_names, drop = FALSE]
    for (jn in names(fi$xlevels)) {
      col <- nd[[jn]]
      tr_lev <- fi$xlevels[[jn]]
      col_chr <- as.character(col)
      if (any(!is.na(col_chr) & !(col_chr %in% tr_lev)))
        stop("ols: out-of-vocabulary factor level(s) in newdata column ",
             jn, ".", call. = FALSE)
      nd[[jn]] <- factor(col, levels = tr_lev)
    }
    xnew <- stats::model.matrix(~ ., data = nd)
    if ("(Intercept)" %in% colnames(xnew))
      xnew <- xnew[, colnames(xnew) != "(Intercept)", drop = FALSE]
    xnew <- xnew[, fi$expanded_names, drop = FALSE]
  } else if (is.data.frame(newdata)) {
    pred_names <- if (has_int) setdiff(cn, "(Intercept)") else cn
    if (!all(pred_names %in% colnames(newdata)))
      stop("ols: newdata is missing columns: ",
           paste(setdiff(pred_names, colnames(newdata)), collapse = ", "),
           call. = FALSE)
    xnew <- as.matrix(newdata[, pred_names, drop = FALSE])
  } else if (is.matrix(newdata) || is.numeric(newdata)) {
    xnew <- as.matrix(newdata)
    n_pred <- if (has_int) ncol(object$X) - 1L else ncol(object$X)
    if (ncol(xnew) != n_pred)
      stop("ols: newdata has ", ncol(xnew),
           " columns; expected ", n_pred, ".", call. = FALSE)
  } else {
    stop("ols: newdata must be a matrix or data frame.", call. = FALSE)
  }
  if (has_int)
    xnew <- cbind(`(Intercept)` = 1, xnew)
  storage.mode(xnew) <- "double"
  xnew
}

# Internal: resolve the coefficient variance-covariance matrix for an
# `ols` fit at the requested robustness level. "none" returns the stored
# classical vcov; HC* recomputes the sandwich via the C++ engine.
#
# @keywords internal
.ols_vcov <- function(object, robust = "none") {
  if (identical(robust, "none")) return(object$vcov)
  type <- switch(robust, HC0 = 1L, HC1 = 2L, HC2 = 3L, HC3 = 4L)
  w_vec <- if (is.null(object$weights)) rep(1, object$n) else
    as.numeric(object$weights)
  v <- ols_vcov_cpp(object$X, w_vec, object$residuals, object$XtXinv,
                    object$hatvalues, object$sigma2, type)
  dimnames(v) <- dimnames(object$vcov)
  v
}

# Internal: per-row observation-level sigma from a stored varmod. Mirrors
# the .ares_make_pi sigma logic but returns the sigma vector directly so
# predict.ols can fold it into a prediction-interval half-width.
#
# @keywords internal
.ols_varmod_sigma <- function(vm, yhat) {
  if (identical(vm$type, "const"))
    return(rep(vm$sigma_hat, length(yhat)))
  # vm$type == "lm": linear MAD model in yhat, floored as in ares.
  raw_mad <- vm$scale * (vm$intercept + vm$slope * yhat)
  floor_val <- if (!is.null(vm$sigma_floor)) vm$sigma_floor else 1e-12
  pmax(raw_mad, floor_val)
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

# ============================================================================
#  predict method for the `logreg` class
# ============================================================================

#' Predictions from a `logreg` fit
#'
#' Returns predictions for new data from a fitted `logreg` model on the
#' link, probability, or class scale, with optional standard errors.
#'
#' @param object An object of class `"logreg"`.
#' @param newdata A numeric matrix or data frame of new predictors.
#'   `NULL` (default) returns predictions for the training data.
#' @param type Prediction scale.
#'   - `"link"`: the linear predictor (log-odds).
#'   - `"response"` (default): the fitted probability.
#'   - `"class"`: a hard 0/1 (or factor) label obtained by thresholding
#'     the fitted probability at `threshold`.
#' @param se.fit If `TRUE`, the returned object carries an `"se.fit"`
#'   attribute holding the per-row standard error. For `type = "link"`
#'   the SE is on the log-odds scale; for `type = "response"` it is the
#'   delta-method SE of the probability. `se.fit` is not available for
#'   `type = "class"`. For bagged fits the SE is the per-row bag
#'   standard deviation of the predicted probability. Default `FALSE`.
#' @param threshold Probability cut-off for `type = "class"`. Default
#'   `0.5`.
#' @param robust Heteroscedasticity-consistent covariance used for the
#'   standard errors.
#'   - `"none"` (default): the classical maximum-likelihood covariance.
#'   - `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`: a sandwich estimator.
#' @param ... Currently unused.
#' @return A numeric vector of predictions on the requested scale (a
#'   factor when `type = "class"` and the fit was built from a factor
#'   response). With `se.fit = TRUE` the result carries an `"se.fit"`
#'   attribute.
#' @examples
#' df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
#' fit <- logreg(y ~ wt + hp, data = df)
#' predict(fit, df[1:5, ], type = "response")
#' predict(fit, df[1:5, ], type = "class")
#' @export
predict.logreg <- function(object, newdata = NULL,
                           type = c("link", "response", "class"),
                           se.fit = FALSE,
                           threshold = 0.5,
                           robust = c("none", "HC0", "HC1", "HC2", "HC3"),
                           ...) {
  type <- match.arg(type)
  robust <- match.arg(robust)
  if (isTRUE(se.fit) && type == "class")
    stop("logreg: se.fit is not available for type = 'class'.",
         call. = FALSE)

  # ---- newdata = NULL fast path -------------------------------------------
  # For a single (non-bagged) fit, the training-data predictions are
  # already stored; with no SE / class request we can return them
  # directly. For a bagged fit, or when SEs are needed, re-route through
  # the design path so the bag mean / vcov can be computed.
  if (is.null(newdata)) {
    if (is.null(object$boot) && !isTRUE(se.fit) && type != "class") {
      if (type == "link") return(object$linear.predictors)
      return(object$fitted.values)              # type == "response".
    }
    xnew <- object$X
  } else {
    xnew <- .logreg_newdata_matrix(object, newdata)
  }
  storage.mode(xnew) <- "double"

  # ---- Linear predictor ---------------------------------------------------
  # Single fit: the central-fit log-odds. Bagged fit: the bag mean is
  # taken on the PROBABILITY scale (averaging probabilities, not
  # log-odds), matching the ares/ols bag-mean convention of averaging the
  # quantity predict() reports.
  eta_central <- drop(xnew %*% object$coefficients)
  boot_probs <- NULL
  if (!is.null(object$boot)) {
    bc <- object$boot$coefficients
    ok <- which(!apply(is.na(bc), 2L, any))
    eta_mat <- cbind(eta_central, xnew %*% bc[, ok, drop = FALSE])
    boot_probs <- stats::plogis(eta_mat)
    prob <- rowMeans(boot_probs)
    eta <- stats::qlogis(pmin(pmax(prob, 1e-10), 1 - 1e-10))
  } else {
    eta <- eta_central
    prob <- stats::plogis(eta)
  }

  # ---- Standard errors ----------------------------------------------------
  se_link <- NULL
  se_resp <- NULL
  if (isTRUE(se.fit)) {
    if (!is.null(object$boot)) {
      # Bagged fit: the per-row SE is the bag standard deviation of the
      # predicted probabilities (plus the central fit).
      se_resp <- apply(boot_probs, 1L, stats::sd)
      se_link <- se_resp / pmax(prob * (1 - prob), 1e-12)
    } else {
      vcov <- .logreg_vcov(object, robust)
      # var(eta) = x' V x, computed rowwise; delta method for the
      # probability scale: d mu / d eta = mu (1 - mu).
      se_link <- sqrt(pmax(rowSums((xnew %*% vcov) * xnew), 0))
      se_resp <- se_link * prob * (1 - prob)
    }
  }

  # ---- Assemble the requested scale --------------------------------------
  if (type == "link") {
    out <- eta
    if (isTRUE(se.fit)) attr(out, "se.fit") <- se_link
    return(out)
  }
  if (type == "response") {
    out <- prob
    if (isTRUE(se.fit)) attr(out, "se.fit") <- se_resp
    return(out)
  }
  # type == "class": threshold the probability. When the fit was built
  # from a factor response, return the original factor levels.
  cls <- as.integer(prob > threshold)
  if (!is.null(object$y.levels))
    cls <- factor(object$y.levels[cls + 1L], levels = object$y.levels)
  cls
}

# Internal: turn `newdata` into a numeric design matrix aligned with the
# training design of a `logreg` fit. Handles the formula path (terms +
# xlevels), the factor-expanded data-frame path, and the plain matrix /
# data-frame path. Mirrors .ols_newdata_matrix.
#
# @keywords internal
.logreg_newdata_matrix <- function(object, newdata) {
  cn <- colnames(object$X)
  has_int <- isTRUE(object$intercept) ||
    (!is.null(cn) && "(Intercept)" %in% cn)
  if (!is.null(object$terms)) {
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, as.data.frame(newdata),
                             xlev = object$xlevels)
    mm <- stats::model.matrix(tt, mf)
    if (!identical(colnames(mm), cn)) {
      missing_cols <- setdiff(cn, colnames(mm))
      if (length(missing_cols))
        stop("logreg: newdata model-matrix expansion is missing columns: ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
      mm <- mm[, cn, drop = FALSE]
    }
    return(mm)
  }
  if (is.data.frame(newdata) && !is.null(object$factor_info)) {
    fi <- object$factor_info
    if (!all(fi$orig_names %in% colnames(newdata)))
      stop("logreg: newdata is missing columns: ",
           paste(setdiff(fi$orig_names, colnames(newdata)),
                 collapse = ", "), call. = FALSE)
    nd <- newdata[, fi$orig_names, drop = FALSE]
    for (jn in names(fi$xlevels)) {
      col <- nd[[jn]]
      tr_lev <- fi$xlevels[[jn]]
      col_chr <- as.character(col)
      if (any(!is.na(col_chr) & !(col_chr %in% tr_lev)))
        stop("logreg: out-of-vocabulary factor level(s) in newdata ",
             "column ", jn, ".", call. = FALSE)
      nd[[jn]] <- factor(col, levels = tr_lev)
    }
    xnew <- stats::model.matrix(~ ., data = nd)
    if ("(Intercept)" %in% colnames(xnew))
      xnew <- xnew[, colnames(xnew) != "(Intercept)", drop = FALSE]
    xnew <- xnew[, fi$expanded_names, drop = FALSE]
  } else if (is.data.frame(newdata)) {
    pred_names <- if (has_int) setdiff(cn, "(Intercept)") else cn
    if (!all(pred_names %in% colnames(newdata)))
      stop("logreg: newdata is missing columns: ",
           paste(setdiff(pred_names, colnames(newdata)), collapse = ", "),
           call. = FALSE)
    xnew <- as.matrix(newdata[, pred_names, drop = FALSE])
  } else if (is.matrix(newdata) || is.numeric(newdata)) {
    xnew <- as.matrix(newdata)
    n_pred <- if (has_int) ncol(object$X) - 1L else ncol(object$X)
    if (ncol(xnew) != n_pred)
      stop("logreg: newdata has ", ncol(xnew),
           " columns; expected ", n_pred, ".", call. = FALSE)
  } else {
    stop("logreg: newdata must be a matrix or data frame.", call. = FALSE)
  }
  if (has_int)
    xnew <- cbind(`(Intercept)` = 1, xnew)
  storage.mode(xnew) <- "double"
  xnew
}

# Internal: resolve the coefficient variance-covariance matrix for a
# `logreg` fit at the requested robustness level. "none" returns the
# stored classical vcov; HC* recomputes the sandwich via the C++ engine.
#
# @keywords internal
.logreg_vcov <- function(object, robust = "none") {
  if (identical(robust, "none")) return(object$vcov)
  type <- switch(robust, HC0 = 1L, HC1 = 2L, HC2 = 3L, HC3 = 4L)
  w_vec <- if (is.null(object$weights)) rep(1, object$n) else
    as.numeric(object$weights)
  v <- logreg_vcov_cpp(object$X, w_vec, object$y, object$fitted.values,
                       object$XtWXinv, object$hatvalues, type)
  dimnames(v) <- dimnames(object$vcov)
  v
}

# ============================================================================
#  predict.bgam -- predictions from a fitted bgam object
# ============================================================================

#' Predictions from a `bgam` fit
#'
#' Returns predictions from a fitted component-wise P-spline boosting model.
#' New B-spline bases are constructed for `newdata` using the stored knot
#' sequences; values outside the training range are extrapolated by the
#' natural boundary behaviour of the de Boor recursion (polynomial
#' extrapolation).  For bagged fits (`n.boot > 0`), the bag mean across all
#' surviving replicates is returned.
#'
#' @param object An object of class `"bgam"`, as returned by [bgam()].
#' @param newdata A numeric matrix or data frame of new predictor values with
#'   the same columns (in any order) as the training data.  `NULL` (default)
#'   returns stored training predictions: `$fitted.values` for non-bagged
#'   fits, bag mean for bagged fits.  Out-of-vocabulary factor levels in
#'   `newdata` cause an informative error.
#' @param type Prediction scale:
#'   \describe{
#'     \item{`"response"`}{(default) Fitted probabilities for
#'       `family = "binomial"`; fitted values (identical to link scale) for
#'       `family = "gaussian"`.}
#'     \item{`"link"`}{The linear predictor (additive predictor on the link
#'       scale).  Identical to `"response"` for gaussian.}
#'     \item{`"terms"`}{An `nrow(newdata)` x `p` numeric matrix of per-
#'       predictor component contributions.  Column `j` holds
#'       `B_new_j %*% beta_j`.  Column names match `object$predictor_names`.
#'       Row sums of all columns plus `object$intercept_value` equal the
#'       link-scale prediction.}
#'   }
#' @param se.fit If `TRUE`, returns a list with elements `$fit` (predictions
#'   as described above) and `$se.fit` (plug-in standard errors summed
#'   across active base-learners).  For predictor j active at observation i,
#'   the component variance is `diag(B_j[i,,drop=FALSE] %*% solve(A_j) %*%
#'   t(B_j[i,,drop=FALSE])) * sigma2` where `A_j = B_j'B_j + lambda_j *
#'   D_j'D_j`.  For binomial, delta-method SE on the response scale is
#'   applied: `se_mu = se_eta * mu_hat * (1 - mu_hat)`.
#'   Default `FALSE`.
#' @param level Unused; reserved for a future prediction-interval interface.
#'   Default `0.95`.
#' @param ... Currently unused.
#'
#' @return
#'   \describe{
#'     \item{`se.fit = FALSE`}{A numeric vector of length `nrow(newdata)`, or
#'       a numeric matrix of dimension `nrow(newdata)` x `p` when
#'       `type = "terms"`.}
#'     \item{`se.fit = TRUE`}{A list with elements `$fit` (as above) and
#'       `$se.fit` (numeric vector of per-observation standard errors).}
#'   }
#'   For bagged fits the bag mean is returned (se.fit is not available for
#'   bagged fits).
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 3), n, 3)
#' y <- sin(x[, 1]) + 0.5 * x[, 2] + rnorm(n, sd = 0.5)
#' fit <- bgam(x, y, mstop = 50, autotune = FALSE)
#' predict(fit, x[1:5, ])
#'
#' # Per-component contributions
#' terms_mat <- predict(fit, x, type = "terms")
#' head(terms_mat)
#'
#' # Plug-in standard errors
#' p_se <- predict(fit, x[1:5, ], se.fit = TRUE)
#' p_se$se.fit
#' @export
predict.bgam <- function(object, newdata = NULL,
                         type   = c("response", "link", "terms"),
                         se.fit = FALSE,
                         level  = 0.95,
                         ...) {
  type <- match.arg(type)

  # ---- NULL newdata: return stored or bag-mean -----------------------------
  if (is.null(newdata)) {
    if (!is.null(object$boot) && length(object$boot$fits) > 0L) {
      # Bagged: route through the training-x path so all boot fits are
      # evaluated at the same n x p training matrix (matches newdata=X path).
      if (!is.null(object$x)) {
        return(Recall(object, newdata = object$x, type = type,
                      se.fit = se.fit, level = level))
      }
      # Fallback if $x was not stored: average fitted.values directly.
      all_fits <- c(list(object), object$boot$fits)
      fv_mat <- do.call(cbind, lapply(all_fits, `[[`, "fitted.values"))
      yhat <- rowMeans(fv_mat, na.rm = TRUE)
      if (isTRUE(se.fit)) {
        se <- .bgam_se_fit(object, newdata_B_list = NULL)
        return(list(fit = yhat, se.fit = se))
      }
      return(yhat)
    }
    # Non-bagged: dispatch on type
    if (identical(type, "terms")) {
      # Build terms matrix from stored training B matrices.
      p <- object$p
      n_obs <- object$n
      terms_mat <- matrix(0, nrow = n_obs, ncol = p,
                          dimnames = list(NULL, object$predictor_names))
      for (j in seq_len(p)) {
        bv <- as.numeric(object$coefficients[[j]])
        if (any(bv != 0)) {
          Bj <- object$base_learners[[j]]$B_train
          if (!is.null(Bj))
            terms_mat[, j] <- as.numeric(Bj %*% bv)
        }
      }
      return(terms_mat)
    }
    yhat <- if (identical(type, "link")) {
      object$linear.predictors
    } else {
      object$fitted.values
    }
    if (isTRUE(se.fit)) {
      se <- .bgam_se_fit(object, newdata_B_list = NULL)
      return(list(fit = yhat, se.fit = se))
    }
    return(yhat)
  }

  # ---- Coerce newdata to matrix -------------------------------------------
  if (is.data.frame(newdata)) {
    if (!is.null(object$terms)) {
      # Formula path: reconstruct model matrix respecting factor levels
      tt <- stats::delete.response(object$terms)
      if (!is.null(object$xlevels)) {
        for (vn in names(object$xlevels)) {
          if (vn %in% names(newdata)) {
            tr_lev <- object$xlevels[[vn]]
            nd_lev <- unique(as.character(newdata[[vn]]))
            bad <- setdiff(nd_lev, c(tr_lev, NA_character_))
            if (length(bad))
              stop("bgam: newdata has factor level(s) not seen in training: ",
                   paste(bad, collapse = ", "), call. = FALSE)
          }
        }
      }
      mf <- stats::model.frame(tt, newdata,
                                xlev = object$xlevels %||% list())
      mm <- stats::model.matrix(tt, mf)
      if ("(Intercept)" %in% colnames(mm))
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      newdata <- mm
    } else {
      # Default path: need columns matching predictor_names
      pn <- object$predictor_names
      if (!all(pn %in% colnames(newdata)))
        stop("bgam: newdata is missing columns: ",
             paste(setdiff(pn, colnames(newdata)), collapse = ", "),
             call. = FALSE)
      newdata <- as.matrix(newdata[, pn, drop = FALSE])
    }
  }
  if (!is.matrix(newdata) || !is.numeric(newdata))
    stop("bgam: newdata must be a numeric matrix or data frame.",
         call. = FALSE)
  storage.mode(newdata) <- "double"

  m_new <- nrow(newdata)
  p <- object$p

  if (ncol(newdata) != p)
    stop("bgam: newdata has ", ncol(newdata), " columns; expected ", p,
         ".", call. = FALSE)

  # ---- Bagged prediction (newdata supplied) --------------------------------
  if (!is.null(object$boot) && !is.null(object$boot$fits)) {
    all_fits <- c(list(object),
                  Filter(Negate(is.null), object$boot$fits))
    preds <- lapply(all_fits, function(f) {
      tryCatch(
        predict.bgam(f, newdata = newdata, type = type, se.fit = FALSE),
        error = function(e) rep(NA_real_, m_new))
    })
    pred_mat <- do.call(cbind, preds)
    return(rowMeans(pred_mat, na.rm = TRUE))
  }

  # ---- Build new B-spline bases for each predictor ------------------------
  B_new_list <- vector("list", p)
  for (j in seq_len(p)) {
    bl <- object$base_learners[[j]]
    xj_new <- newdata[, j]
    if (length(bl$knots) == 2L && bl$K == 1L) {
      # Linear base-learner (near-constant or unpenalized)
      B_new_list[[j]] <- matrix(xj_new - mean(xj_new), ncol = 1L)
    } else {
      res <- bgam_bspline_basis_cpp(xj_new, bl$knots, bl$degree)
      B_new_list[[j]] <- res$B
    }
  }

  # ---- Link-scale prediction ----------------------------------------------
  eta_new <- bgam_predict_cpp(B_new_list, object$coefficients,
                               object$intercept_value)

  # ---- Terms: per-component matrix ----------------------------------------
  if (identical(type, "terms")) {
    terms_mat <- matrix(0, nrow = m_new, ncol = p,
                        dimnames = list(NULL, object$predictor_names))
    for (j in seq_len(p)) {
      bv <- as.numeric(object$coefficients[[j]])
      if (any(bv != 0))
        terms_mat[, j] <- as.numeric(B_new_list[[j]] %*% bv)
    }
    return(terms_mat)
  }

  # ---- Response scale ------------------------------------------------------
  if (identical(type, "link")) {
    yhat <- as.numeric(eta_new)
  } else {
    yhat <- if (identical(object$family, "binomial")) {
      as.numeric(1 / (1 + exp(-eta_new)))
    } else {
      as.numeric(eta_new)
    }
  }

  # ---- se.fit (plug-in) ---------------------------------------------------
  if (isTRUE(se.fit)) {
    se <- .bgam_se_fit(object, newdata_B_list = B_new_list)
    if (identical(object$family, "binomial") && !identical(type, "link")) {
      mu_hat <- 1 / (1 + exp(-as.numeric(eta_new)))
      se <- se * mu_hat * (1 - mu_hat)
    }
    return(list(fit = yhat, se.fit = se))
  }

  yhat
}

# Internal: compute plug-in SE for bgam predictions.
# Formula: SE_i = sqrt( sum_{j active} diag(B_j[i,] A_j^{-1} B_j[i,]') * sigma2 )
# where A_j = B_j'B_j + lambda_j P_j and L_j = chol(A_j, "lower").
# Equiv: sum rowSums( (B_j %*% t(solve(L_j)))^2 ) * sigma2
.bgam_se_fit <- function(object, newdata_B_list = NULL) {
  p <- object$p
  sigma2 <- object$sigma2
  if (is.null(sigma2)) sigma2 <- 1  # binomial: SE on link scale

  # Determine per-predictor B matrices to use.
  if (is.null(newdata_B_list)) {
    # Use stored training B matrices for per-observation SE on training data.
    B_use <- lapply(object$base_learners, `[[`, "B_train")
  } else {
    B_use <- newdata_B_list
  }
  n_out <- nrow(B_use[[1L]])

  var_total <- numeric(n_out)
  for (j in seq_len(p)) {
    bv <- as.numeric(object$coefficients[[j]])
    if (!any(bv != 0)) next
    bl <- object$base_learners[[j]]
    L_j <- bl$chol  # K_j x K_j lower Cholesky of A_j = B_j'B_j + lambda_j P_j
    Bj  <- B_use[[j]]
    # solve(t(L_j), I) = L_j^{-T}: backsolve on upper-triangular t(L_j)
    Li_inv <- backsolve(t(L_j), diag(ncol(L_j)))  # upper.tri=TRUE (default)
    # (B_j L_j^{-T}) so rowSums(.^2) = diag(B_j A_j^{-1} B_j')
    Bj_Linv <- Bj %*% t(Li_inv)
    var_j   <- rowSums(Bj_Linv^2) * sigma2
    var_total <- var_total + var_j
  }
  sqrt(pmax(var_total, 0))
}

# R-level lower-triangular Cholesky solve helper.
.bgam_chol_solve_r <- function(L, b) {
  y <- forwardsolve(L, b)
  backsolve(t(L), y)
}
