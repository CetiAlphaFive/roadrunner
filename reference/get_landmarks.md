# Extract landmark coordinates from a Nystrom krls fit

Returns the m x d matrix of landmark coordinates used by a fit built
with `approx = "nystrom"`. Internally landmarks live in standardized
X-space; the default `scale = "original"` undoes the standardization
using the training X's centers and SDs.

## Usage

``` r
get_landmarks(fit, scale = c("original", "standardized"))
```

## Arguments

- fit:

  A `krls_rr` fit built with `approx = "nystrom"`.

- scale:

  Either `"original"` (default) or `"standardized"`.

## Value

Numeric matrix of dimension `nystrom_m` x `ncol(X_train)`.
