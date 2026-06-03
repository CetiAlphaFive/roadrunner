# Diagnostic plots for a `meep` cross-fitted ensemble

Visualises out-of-fold (OOF) ensemble performance per nuisance. For each
selected nuisance the panel layout is keyed to that nuisance's family:

- **binomial** – one ROC panel overlaying each present learner's OOF ROC
  curve (thin, coloured) plus the combined "stack" curve (thick, black),
  with a 45-degree reference line and a legend reporting each curve's
  AUC.

- **gaussian** – two panels: an OOF \\R^2\\ bar chart (one bar per
  present learner plus a "stack" bar) and an observed-vs-OOF-prediction
  scatter for the stack with a 45-degree line and the stack \\R^2\\
  annotated.

## Usage

``` r
# S3 method for class 'meep'
plot(x, which = c("outcome", "treatment"), ...)
```

## Arguments

- x:

  An object of class `"meep"`, as returned by
  [`meep()`](https://cetialphafive.github.io/roadrunner/reference/meep.md).

- which:

  Character vector of nuisances to draw; any subset of
  `c("outcome", "treatment", "mu0", "mu1")`. Nuisances not present on
  the object are silently dropped. Defaults to
  `c("outcome", "treatment")`.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.

## Details

Panels are laid out on a near-square grid sized to the total panel count
(binomial contributes 1 panel, gaussian 2). The function works in
headless (non-interactive) mode without error; degenerate nuisances
(single-class response, zero-variance outcome) are skipped gracefully
rather than erroring.

## See also

[`meep()`](https://cetialphafive.github.io/roadrunner/reference/meep.md)
