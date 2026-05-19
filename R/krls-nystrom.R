# Internal helpers for krls(..., approx = "nystrom").
# Phase 2 (v0.0.0.9043). See inst/specs/2026-05-19-krls-nystrom-design.md.
#
# Nystrom replaces the full n x n kernel by a low-rank approximation
# anchored at m << n landmarks. The ridge problem then lives in
# m-dimensional space:
#
#   Phi      = C %*% U %*% diag(D_reg^{-1/2})       (n x m feature map)
#   beta_hat = (Phi' Phi + lambda I)^{-1} Phi' y
#   f_hat(x) = K(x, Z) %*% alpha,   alpha = U %*% (D_reg^{-1/2} * beta)
#
# where C = K(X, Z), W = K(Z, Z), and (U, D) is the eigendecomposition
# of W with relative-ridge stabilization D_reg = pmax(D, eps * max(D)).

.select_landmarks_random <- function(n, m) {
  sort.int(sample.int(n, m))
}

# Run `expr` with the RNG temporarily seeded at `seed`, restoring the
# caller's .Random.seed (or its absence) on exit. Used so that
# passing `landmark_seed` to krls() does NOT perturb the user's
# downstream RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  expr
}

.validate_nystrom_m <- function(nystrom_m, n) {
  if (!is.numeric(nystrom_m) || length(nystrom_m) != 1L ||
      is.na(nystrom_m) || !is.finite(nystrom_m) ||
      nystrom_m != as.integer(nystrom_m)) {
    stop("nystrom_m must be a finite integer")
  }
  nystrom_m <- as.integer(nystrom_m)
  if (nystrom_m < 1L || nystrom_m > n) {
    stop("nystrom_m must satisfy 1 <= nystrom_m <= nrow(X)")
  }
  nystrom_m
}

.validate_nystrom_eps <- function(nystrom_eps) {
  if (!is.numeric(nystrom_eps) || length(nystrom_eps) != 1L ||
      is.na(nystrom_eps) || !is.finite(nystrom_eps) ||
      nystrom_eps <= 0) {
    stop("nystrom_eps must be a finite positive scalar")
  }
  nystrom_eps
}

# X_std is the already-standardized training matrix; X_centers / X_scales
# are the column means / SDs of the original-scale training X.
# Index- and NULL-form landmarks resolve via X_std rows (already standardized).
# Matrix-form landmarks are taken as user-supplied original-scale
# coordinates and standardized with the same centers/scales.
.resolve_landmarks <- function(landmarks, landmark_method, nystrom_m,
                               X_std, X_centers, X_scales,
                               landmark_seed = NULL) {
  n <- nrow(X_std)
  d <- ncol(X_std)

  if (is.null(landmarks)) {
    if (is.null(nystrom_m)) {
      nystrom_m <- ceiling(sqrt(n) * 2)
    }
    nystrom_m <- .validate_nystrom_m(nystrom_m, n)
    return(.with_seed(landmark_seed, {
      if (identical(landmark_method, "kmeans")) {
        if (nystrom_m == n) {
          list(indices = seq_len(n),
               matrix  = X_std,
               method_used = "kmeans")
        } else {
          km <- stats::kmeans(X_std, centers = nystrom_m,
                              nstart = 10L, iter.max = 50L)
          Z  <- unname(km$centers)
          list(indices = NULL, matrix = Z, method_used = "kmeans")
        }
      } else {
        idx <- .select_landmarks_random(n, nystrom_m)
        list(indices = idx,
             matrix  = X_std[idx, , drop = FALSE],
             method_used = "random")
      }
    }))
  }

  if (is.numeric(landmarks) && is.null(dim(landmarks))) {
    if (anyNA(landmarks) || any(!is.finite(landmarks)) ||
        any(landmarks != as.integer(landmarks))) {
      stop("landmarks (as indices) must be integer-valued")
    }
    idx <- as.integer(landmarks)
    if (anyNA(idx) || any(idx < 1L) || any(idx > n) || anyDuplicated(idx)) {
      stop("landmarks (as indices) must be unique integers in 1:nrow(X)")
    }
    return(list(indices = idx,
                matrix  = X_std[idx, , drop = FALSE],
                method_used = "user_indices"))
  }

  if (is.matrix(landmarks) || is.data.frame(landmarks)) {
    Z <- as.matrix(landmarks)
    if (!is.numeric(Z))       stop("landmarks matrix must be numeric")
    if (ncol(Z) != d)         stop("ncol(landmarks) must equal ncol(X)")
    if (anyNA(Z) || any(!is.finite(Z)))
      stop("landmarks matrix must contain only finite, non-NA values")
    if (anyDuplicated(Z))
      stop("landmarks matrix must not contain duplicate rows ",
           "(W would be rank-deficient)")
    Z_std <- scale(Z, center = X_centers, scale = X_scales)
    attr(Z_std, "scaled:center") <- NULL
    attr(Z_std, "scaled:scale")  <- NULL
    return(list(indices = NULL, matrix = Z_std, method_used = "user_matrix"))
  }

  stop("`landmarks` must be NULL, an integer vector of row indices, ",
       "or an m x d numeric matrix (in original X-scale)")
}
