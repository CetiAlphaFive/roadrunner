# Package index

## MARS: ares()

Fast MARS with GCV / CV pruning, autotune, bagging, weights, families,
and prediction intervals.

- [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  : Fast Multivariate Adaptive Regression Splines

- [`predict(`*`<ares>`*`)`](https://cetialphafive.github.io/roadrunner/reference/predict.ares.md)
  :

  Predictions from an `ares` fit

- [`print(`*`<ares>`*`)`](https://cetialphafive.github.io/roadrunner/reference/print.ares.md)
  :

  Print method for `ares` fits

- [`summary(`*`<ares>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.ares.md)
  [`print(`*`<summary.ares>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.ares.md)
  :

  Summary method for `ares` fits

- [`plot(`*`<ares>`*`)`](https://cetialphafive.github.io/roadrunner/reference/plot.ares.md)
  :

  Diagnostic plots for an `ares` fit

## KRLS: krls()

Kernel Regularized Least Squares with closed-form LOO lambda selection
and marginal effects.

- [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  [`predict(`*`<krls_rr>`*`)`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  : Kernel Regularized Least Squares

- [`plot(`*`<krls_rr>`*`)`](https://cetialphafive.github.io/roadrunner/reference/plot.krls_rr.md)
  : Diagnostic plot for KRLS fits

- [`slim_krls()`](https://cetialphafive.github.io/roadrunner/reference/slim_krls.md)
  : Strip heavy fields off a KRLS fit

- [`unslim_krls()`](https://cetialphafive.github.io/roadrunner/reference/unslim_krls.md)
  :

  Undo
  [`slim_krls()`](https://cetialphafive.github.io/roadrunner/reference/slim_krls.md)
  (placeholder)

- [`get_landmarks()`](https://cetialphafive.github.io/roadrunner/reference/get_landmarks.md)
  : Extract landmark coordinates from a Nystrom krls fit

## Causal ensemble: meep()

Cross-fitted ensemble of ares() and krls() for DML and causal-forest
nuisance estimation.

- [`meep()`](https://cetialphafive.github.io/roadrunner/reference/meep.md)
  :

  Cross-fitted causal ensemble of
  [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  and
  [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)

- [`predict(`*`<meep>`*`)`](https://cetialphafive.github.io/roadrunner/reference/predict.meep.md)
  :

  Predict from a fitted `meep` object

## Package

- [`roadrunner`](https://cetialphafive.github.io/roadrunner/reference/roadrunner-package.md)
  [`roadrunner-package`](https://cetialphafive.github.io/roadrunner/reference/roadrunner-package.md)
  : roadrunner: Fast, Low-Dependency Machine Learning Algorithms
