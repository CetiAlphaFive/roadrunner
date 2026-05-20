# Predict from a fitted `meep` object

Out-of-fold predictions exist only for the training rows (and are stored
on the fitted object). For genuinely new data, `predict.meep()` refits
each base learner on the **full** training sample – using the frozen
hyperparameters when `tune = "once"` – then combines the refits with the
stored ensemble weights. Full-data refits are cached lazily in an
environment field on the object and are not serialized.

## Usage

``` r
# S3 method for class 'meep'
predict(
  object,
  newdata,
  nuisance = c("outcome", "treatment", "mu0", "mu1"),
  ...
)
```

## Arguments

- object:

  A fitted `"meep"` object.

- newdata:

  A matrix or data.frame of new covariates.

- nuisance:

  Which nuisance to predict: `"outcome"`, `"treatment"`, `"mu0"`, or
  `"mu1"`.

- ...:

  Unused.

## Value

A numeric vector of length `nrow(newdata)`.

## Details

Predictions on new data can extrapolate; the cross-fitting honesty
guarantee applies only to the training-row OOF predictions.
