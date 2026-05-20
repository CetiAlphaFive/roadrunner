# Strip heavy fields off a KRLS fit

Removes large `n x n` (or `n x m`) intermediates from a fitted `krls_rr`
object so it serialises (`saveRDS`) much smaller. Two modes:

## Usage

``` r
slim_krls(fit, keep_predict = TRUE)
```

## Arguments

- fit:

  A `krls_rr` fit.

- keep_predict:

  Logical. If `TRUE` (default), keep enough state for
  [`predict()`](https://rdrr.io/r/stats/predict.html) to run.

## Value

A `krls_rr` object with heavy fields removed and `$slimmed = TRUE`.

## Details

- `keep_predict = TRUE` (default): keeps everything
  [`predict.krls_rr()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  needs (coeffs, sigma / sigma_vec, lambda, the standardisation stats,
  factor bookkeeping, and Nystrom landmarks where applicable). Drops the
  `K`, eigenvector matrices `V` / `Vsq`, projected response `Vty`, and
  `vcov.c`.

- `keep_predict = FALSE`: keeps only inspection-grade summary fields
  (`coeffs`, `R2`, `Looe`, `avgderivatives`, `var.avgderivatives`,
  `sigma`, `sigma_vec`, `lambda`, `Neffective`).
  [`predict()`](https://rdrr.io/r/stats/predict.html) will refuse to
  run.

The returned object carries the flag `$slimmed = TRUE`.
[`unslim_krls()`](https://cetialphafive.github.io/roadrunner/reference/unslim_krls.md)
is a no-op (heavy intermediates can only be rebuilt by refitting).
