# Undo `slim_krls()` (placeholder)

Heavy fields (`K`, `V`, `Vsq`, `Vty`, `vcov.c`, `vcov.fitted`) cannot be
reconstructed without `X` and `y` and the original hyperparameters. This
function therefore re-fits is *not* attempted; it simply returns the
input unchanged and (optionally) prints a hint to refit. Reserved for a
future release where re-construction from cached state may become
available.

## Usage

``` r
unslim_krls(fit)
```

## Arguments

- fit:

  A `krls_rr` fit previously passed through
  [`slim_krls()`](https://cetialphafive.github.io/roadrunner/reference/slim_krls.md).

## Value

`fit`, unchanged. Refit
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
with the original call to rebuild the dropped fields.
