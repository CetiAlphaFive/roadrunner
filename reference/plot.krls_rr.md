# Diagnostic plot for KRLS fits

Four-panel diagnostic display (residuals vs fitted, Q-Q on standardised
residuals, scale-location, residuals vs leverage with Cook contours) in
a 2x2 grid by default. Panels 4 (Cook's distance) and 6 (Cook vs
leverage) are available via `which = 1:6`. Modelled on
[`stats::plot.lm()`](https://rdrr.io/r/stats/plot.lm.html) and
[`plot.ares()`](https://cetialphafive.github.io/roadrunner/reference/plot.ares.md).

## Usage

``` r
# S3 method for class 'krls_rr'
plot(
  x,
  which = c(1L, 2L, 3L, 5L),
  caption = list("Residuals vs Fitted", "Normal Q-Q", "Scale-Location",
    "Cook's distance", "Residuals vs Leverage", expression("Cook's dist vs Leverage  " *
    h[ii]/(1 - h[ii]))),
  panel = if (add.smooth) graphics::panel.smooth else graphics::points,
  sub.caption = NULL,
  main = "",
  ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
  ...,
  id.n = 3L,
  labels.id = NULL,
  cex.id = 0.75,
  qqline = TRUE,
  cook.levels = c(0.5, 1),
  add.smooth = getOption("add.smooth", TRUE),
  label.pos = c(4, 2)
)
```

## Arguments

- x:

  A fitted `krls_rr` object.

- which:

  Which panels to draw; integer subset of `1:6`. Defaults to
  `c(1L, 2L, 3L, 5L)`.

- caption:

  Panel titles (length 6).

- panel:

  Per-panel display function (defaults to LOESS smoother).

- sub.caption:

  Optional sub-caption shown across the figure.

- main:

  Per-panel title.

- ask:

  Activate interactive panel paging.

- ...:

  Passed to underlying `plot` calls.

- id.n:

  Number of extreme points to label per panel.

- labels.id:

  Optional vector of point labels.

- cex.id:

  Label cex.

- qqline:

  Draw the Q-Q line.

- cook.levels:

  Cook's-distance contour levels.

- add.smooth:

  Add LOESS smoother overlay.

- label.pos:

  Label position codes.

## Value

Invisibly returns `x`.

## Details

Leverage uses the KRLS hat matrix \$\$H = K (K + \lambda I)^{-1} = V
\mathrm{diag}(d/(d+\lambda)) V'\$\$ so
`h_i = sum_k V[i,k]^2 * d_k / (d_k + lambda)`. The effective degrees of
freedom is `sum(d/(d+lambda))`. Residual dispersion is
`sqrt(RSS / (n - effdf))` (or `out$varmod$sigma_hat` when present).
