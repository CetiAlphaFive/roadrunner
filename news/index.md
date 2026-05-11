# Changelog

## roadrunner 0.0.0.9028 (development)

### Package rename: `ares` -\> `roadrunner`

The package is rebranded as `roadrunner` – a collection of fast,
low-dependency machine learning algorithm implementations. The MARS
fitter stays exposed as
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md).
Existing
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
call sites, the `"ares"` S3 class,
[`predict.ares()`](https://cetialphafive.github.io/roadrunner/reference/predict.ares.md),
[`print.ares()`](https://cetialphafive.github.io/roadrunner/reference/print.ares.md),
[`summary.ares()`](https://cetialphafive.github.io/roadrunner/reference/summary.ares.md),
and
[`plot.ares()`](https://cetialphafive.github.io/roadrunner/reference/plot.ares.md)
are unchanged. Only the package itself is renamed:

- [`library(ares)`](https://rdrr.io/r/base/library.html) -\>
  [`library(roadrunner)`](https://cetialphafive.github.io/roadrunner/)
- `ares::ares(...)` -\> `roadrunner::ares(...)`
- `Package: ares` in DESCRIPTION -\> `Package: roadrunner`

No user-visible behavior in the MARS fitter changes at this version.
