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

## Penalized LDA: plda()

Penalized Fisher’s linear discriminant (Witten & Tibshirani 2011) with
L1 and fused-lasso penalties, multi-class, and CV autotune.

- [`plda()`](https://cetialphafive.github.io/roadrunner/reference/plda.md)
  : Penalized linear discriminant analysis

- [`predict(`*`<plda>`*`)`](https://cetialphafive.github.io/roadrunner/reference/predict.plda.md)
  : Predict from a penalized LDA fit

- [`print(`*`<plda>`*`)`](https://cetialphafive.github.io/roadrunner/reference/print.plda.md)
  :

  Print method for `plda` fits

- [`summary(`*`<plda>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.plda.md)
  [`print(`*`<summary.plda>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.plda.md)
  :

  Summary method for `plda` fits

- [`plot(`*`<plda>`*`)`](https://cetialphafive.github.io/roadrunner/reference/plot.plda.md)
  :

  Projection plot for `plda` fits

## Linear models: ols()

Ordinary and weighted least squares with classical and HC robust
standard errors, bagging, and prediction intervals.

- [`ols()`](https://cetialphafive.github.io/roadrunner/reference/ols.md)
  : Ordinary and weighted least squares

- [`predict(`*`<ols>`*`)`](https://cetialphafive.github.io/roadrunner/reference/predict.ols.md)
  :

  Predictions from an `ols` fit

- [`print(`*`<ols>`*`)`](https://cetialphafive.github.io/roadrunner/reference/print.ols.md)
  :

  Print method for `ols` fits

- [`summary(`*`<ols>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.ols.md)
  [`print(`*`<summary.ols>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.ols.md)
  :

  Summary method for `ols` fits

- [`plot(`*`<ols>`*`)`](https://cetialphafive.github.io/roadrunner/reference/plot.ols.md)
  :

  Diagnostic plots for an `ols` fit

## Logistic regression: logreg()

Binary logistic regression by IRLS with classical and HC robust standard
errors and bagging.

- [`logreg()`](https://cetialphafive.github.io/roadrunner/reference/logreg.md)
  : Binary logistic regression

- [`predict(`*`<logreg>`*`)`](https://cetialphafive.github.io/roadrunner/reference/predict.logreg.md)
  :

  Predictions from a `logreg` fit

- [`print(`*`<logreg>`*`)`](https://cetialphafive.github.io/roadrunner/reference/print.logreg.md)
  :

  Print method for `logreg` fits

- [`summary(`*`<logreg>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.logreg.md)
  [`print(`*`<summary.logreg>`*`)`](https://cetialphafive.github.io/roadrunner/reference/summary.logreg.md)
  :

  Summary method for `logreg` fits

- [`plot(`*`<logreg>`*`)`](https://cetialphafive.github.io/roadrunner/reference/plot.logreg.md)
  :

  Diagnostic plots for a `logreg` fit

## Package

- [`roadrunner`](https://cetialphafive.github.io/roadrunner/reference/roadrunner-package.md)
  [`roadrunner-package`](https://cetialphafive.github.io/roadrunner/reference/roadrunner-package.md)
  : roadrunner: Fast, Low-Dependency Machine Learning Algorithms
