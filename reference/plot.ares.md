# Diagnostic plots for an `ares` fit

Six diagnostic panels modelled on
[`stats::plot.lm()`](https://rdrr.io/r/stats/plot.lm.html): residuals vs
fitted, normal Q-Q of standardized residuals, scale-location, Cook's
distance, residuals vs leverage (with Cook's-distance contours), and
Cook's distance vs leverage. The default `which` draws panels 1, 2, 3, 5
in a 2x2 grid when the device is on a single-panel layout; this matches
the `plot.lm()` default. For non-gaussian families (`binomial`,
`poisson`, `gamma`), deviance residuals are used in panel 1 and IRLS
working weights are used to compute leverages and standardized residuals
(the GLM analogue of the gaussian hat matrix).

## Usage

``` r
# S3 method for class 'ares'
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

  A fitted object of class `"ares"`.

- which:

  Integer subset of `1:6` selecting which panels to draw. Default
  `c(1, 2, 3, 5)` (the `plot.lm()` default).

- caption:

  Length-6 list of panel captions; defaults match `plot.lm()`.

- panel:

  Panel function for scatter panels; defaults to
  [`graphics::panel.smooth()`](https://rdrr.io/r/graphics/panel.smooth.html).

- sub.caption:

  Optional caption shown on the outer margin (one per page). Defaults to
  a deparsed `x$call`.

- main:

  Title for individual panels (typically `""`).

- ask:

  If `TRUE`, prompt before each page when multiple panels are drawn on a
  single-panel device.

- ...:

  Further graphical parameters passed through to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) calls.

- id.n:

  Number of extreme points to label per panel (default `3`; set `0` to
  suppress).

- labels.id:

  Character vector of row labels; defaults to `1:n`.

- cex.id:

  `cex` for point labels.

- qqline:

  Logical; if `TRUE`, add a Q-Q reference line on panel 2.

- cook.levels:

  Cook's-distance contours drawn on panel 5.

- add.smooth:

  Logical; passed to `panel`.

- label.pos:

  `pos` argument forwarded to
  [`graphics::text()`](https://rdrr.io/r/graphics/text.html).

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
  plot(fit)
} # }
```
