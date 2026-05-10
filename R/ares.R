# ares — main fitting function

#' Fast Multivariate Adaptive Regression Splines
#'
#' Fits a MARS model (Friedman 1991) using a fast least-squares forward pass
#' with a parallelized knot search and a backward subset-selection that
#' minimizes the GCV criterion. The implementation aims for numerical parity
#' with [earth::earth()] on the gaussian-only core while taking advantage of
#' multi-core CPUs to reduce wall-clock fitting time.
#'
#' Two interfaces are provided. `ares.default()` accepts a numeric predictor
#' matrix `x` and a numeric response vector `y`. `ares.formula()` accepts a
#' model formula and a data frame.
#'
#' @param x A numeric matrix or data frame of predictors, OR a model formula
#'   when calling the formula method.
#' @param y A numeric response vector. Length must equal `nrow(x)`.
#' @param data Used by the formula method only. A data frame.
#' @param degree Maximum interaction degree. Default 1.
#' @param nk Maximum number of basis terms in the forward pass. Default
#'   `min(200, max(20, 2 * ncol(x))) + 1`.
#' @param penalty GCV penalty. Default `if (degree > 1) 3 else 2`.
#' @param thresh Forward-pass relative-RSS early-stop threshold. Default 0.001.
#' @param minspan Minimum knot span. 0 (default) selects an automatic value.
#' @param endspan Knot offset from the data ends. 0 (default) selects automatic.
#' @param adjust.endspan Multiplier applied to `endspan` when the candidate
#'   hinge would deepen an existing interaction (parent term has degree
#'   >= 1). Default `1` (no adjustment). Earth's analogous default is `2`;
#'   ares's heuristic is empirically mixed on the inst/sims grid (helps a
#'   few cells, hurts others), so the default is conservative. Pass `2`
#'   to enable.
#' @param auto.linpreds If `TRUE`, and the best forward-pass knot for a
#'   candidate hinge sits at the boundary of its variable's eligible range,
#'   substitutes a linear term (`dirs = 2`) for the hinge pair. **Default
#'   `FALSE`.** Earth's analogous default is `TRUE`, but earth applies a
#'   stricter linearity test that ares's "boundary knot" heuristic only
#'   approximates; current ares implementation regresses parity on most
#'   benchmark cells. Treat as experimental.
#' @param nprune Maximum number of terms after backward pruning (default `nk`).
#' @param pmethod Pruning method: `"backward"` (default) or `"none"`.
#' @param trace Trace level. 0 = silent (default), 1 = forward-pass progress.
#' @param nthreads Number of threads. 0 (default) selects
#'   `RcppParallel::defaultNumThreads()`. CRAN-distributed examples and the
#'   bundled vignette cap this at 2.
#' @param ... Additional arguments. Currently ignored.
#' @return An object of class `"ares"` — a list with components
#'   `coefficients`, `bx`, `dirs`, `cuts`, `selected.terms`, `rss`, `gcv`,
#'   `rss.per.subset`, `gcv.per.subset`, `fitted.values`, `residuals`, `namesx`,
#'   `call`, plus echoed control parameters.
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' *Annals of Statistics* 19(1):1-67.
#' @examples
#' fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
#' print(fit)
#' p <- predict(fit, as.matrix(mtcars[, -1]))
#' @export
ares <- function(x, ...) UseMethod("ares")

#' @rdname ares
#' @export
ares.formula <- function(x, data = NULL, ..., y = NULL) {
  formula <- x
  if (is.null(data)) data <- environment(formula)
  cl <- match.call()
  mf <- stats::model.frame(formula, data = data)
  yv <- stats::model.response(mf)
  if (is.null(yv)) stop("ares: response variable is missing from formula/data.")
  mm <- stats::model.matrix(formula, mf)
  # drop intercept column
  has_int <- "(Intercept)" %in% colnames(mm)
  if (has_int) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  out <- ares.default(x = mm, y = as.numeric(yv), ...)
  out$call <- cl
  out$terms <- stats::terms(formula, data = data)
  out
}

#' @rdname ares
#' @export
ares.default <- function(x, y, degree = 1L, nk = NULL, penalty = NULL,
                         thresh = 0.001, minspan = 0L, endspan = 0L,
                         adjust.endspan = 1L, auto.linpreds = FALSE,
                         fast.k = 10L, fast.beta = 1.0,
                         nprune = NULL, pmethod = c("backward", "none"),
                         trace = 0L, nthreads = 0L, ...) {
  cl <- match.call()
  pmethod <- match.arg(pmethod)
  pmethod_int <- if (pmethod == "backward") 0L else 1L

  # ---- Input validation ----
  if (is.data.frame(x)) {
    nm <- names(x)
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    colnames(x) <- nm
  }
  if (!is.matrix(x) || !is.numeric(x))
    stop("ares: x must be a numeric matrix or data frame.")
  if (any(!is.finite(x)))
    stop("ares: x contains NA, NaN, or non-finite values.")
  if (!is.numeric(y))
    stop("ares: y must be numeric.")
  if (any(!is.finite(y)))
    stop("ares: y contains NA, NaN, or non-finite values.")
  if (length(y) != nrow(x))
    stop("ares: length(y) (", length(y), ") must equal nrow(x) (", nrow(x), ").")
  if (!is.numeric(degree) || length(degree) != 1L || degree < 1)
    stop("ares: degree must be a single integer >= 1.")
  degree <- as.integer(degree)
  if (ncol(x) < 1L)
    stop("ares: x must have at least one column.")

  # ---- Drop constant columns with warning ----
  col_const <- vapply(seq_len(ncol(x)), function(j) {
    z <- x[, j]
    diff(range(z)) < .Machine$double.eps * 16
  }, logical(1))
  dropped_names <- character(0)
  if (any(col_const)) {
    dropped_names <- colnames(x)[col_const]
    if (is.null(dropped_names)) dropped_names <- paste0("V", which(col_const))
    warning("ares: dropped ", sum(col_const), " constant column(s): ",
            paste(dropped_names, collapse = ", "), call. = FALSE)
    x <- x[, !col_const, drop = FALSE]
    if (ncol(x) < 1L)
      stop("ares: all predictor columns are constant after dropping.")
  }

  # ---- Defaults ----
  p <- ncol(x)
  if (is.null(nk)) nk <- min(200L, max(20L, 2L * p)) + 1L
  nk <- as.integer(nk)
  if (nk < 3L) {
    warning("ares: nk coerced to 3 (minimum).", call. = FALSE)
    nk <- 3L
  }
  if (is.null(penalty)) penalty <- if (degree > 1L) 3 else 2
  penalty <- as.numeric(penalty)
  if (is.null(nprune)) nprune <- nk
  nprune <- as.integer(nprune)
  trace <- as.integer(trace)
  minspan <- as.integer(minspan)
  endspan <- as.integer(endspan)
  adjust_endspan <- as.integer(adjust.endspan)
  if (adjust_endspan < 1L) {
    warning("ares: adjust.endspan coerced to 1 (minimum).", call. = FALSE)
    adjust_endspan <- 1L
  }
  auto_linpreds <- as.integer(isTRUE(auto.linpreds))
  fast_k <- as.integer(fast.k)
  if (is.na(fast_k) || fast_k < 0L) fast_k <- 0L  # 0 = unlimited (no caching)
  fast_beta <- as.numeric(fast.beta)
  if (is.na(fast_beta) || fast_beta < 0) fast_beta <- 0

  nthreads <- as.integer(nthreads)
  nthreads_eff <- if (nthreads <= 0L) RcppParallel::defaultNumThreads() else nthreads
  RcppParallel::setThreadOptions(numThreads = nthreads_eff)

  # ---- Capture column names ----
  namesx <- colnames(x)
  if (is.null(namesx)) namesx <- paste0("V", seq_len(p))

  # ---- Coerce x to plain numeric matrix ----
  storage.mode(x) <- "double"

  # ---- Call C++ engine ----
  out <- mars_fit_cpp(x, as.numeric(y), degree, nk, penalty, thresh,
                      minspan, endspan, adjust_endspan, auto_linpreds,
                      fast_k, fast_beta,
                      nprune, pmethod_int, trace, nthreads_eff)

  # ---- Post-process ----
  rownames(out$dirs) <- .term_labels(out$dirs, out$cuts, namesx)
  rownames(out$cuts) <- rownames(out$dirs)
  colnames(out$dirs) <- namesx
  colnames(out$cuts) <- namesx
  colnames(out$bx) <- rownames(out$dirs)[out$selected.terms]
  names(out$coefficients) <- colnames(out$bx)
  out$fitted.values <- drop(out$bx %*% out$coefficients)
  out$residuals <- as.numeric(y) - out$fitted.values
  out$namesx <- namesx
  out$call <- cl
  out$pmethod <- pmethod
  out$dropped <- dropped_names
  class(out) <- c("ares")
  out
}

# Internal: build term labels matching earth's "h(x-cut)" format
.term_labels <- function(dirs, cuts, namesx) {
  M <- nrow(dirs)
  p <- ncol(dirs)
  labs <- character(M)
  labs[1] <- "(Intercept)"
  for (t in seq_len(M)[-1L]) {
    parts <- character(0)
    for (j in seq_len(p)) {
      d <- dirs[t, j]
      if (d == 0) next
      cv <- cuts[t, j]
      if (d == 2) {
        parts <- c(parts, namesx[j])
      } else if (d == 1) {
        parts <- c(parts, sprintf("h(%s-%g)", namesx[j], cv))
      } else if (d == -1) {
        parts <- c(parts, sprintf("h(%g-%s)", cv, namesx[j]))
      }
    }
    labs[t] <- if (length(parts)) paste(parts, collapse = " * ") else paste0("term", t)
  }
  labs
}
