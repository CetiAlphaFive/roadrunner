# roadrunner: Fast, Low-Dependency Machine Learning Algorithms

`roadrunner` is a collection of fast, low-dependency implementations of
classical machine learning algorithms with thin, base-R-style
interfaces. The two algorithms shipped today are:

## Details

- [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  – Multivariate Adaptive Regression Splines (MARS) with GCV /
  cross-validated pruning, hyperparameter autotune, bagging, weights,
  binomial / poisson / gamma families, and prediction intervals. Built
  on Friedman's (1991) fast least-squares forward selection with a
  parallelised knot search.

- [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  – Kernel Regularized Least Squares (Hainmueller and Hazlett 2014) with
  closed-form leave-one-out lambda selection, per-observation marginal
  effects, and binary first-difference handling.

Both engines are C++ on top of `Rcpp`, `RcppArmadillo`, and
`RcppParallel`. Fits are deterministic across thread counts at a fixed
seed.

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics*, 19(1):1-67.

Hainmueller, J. and Hazlett, C. (2014). Kernel Regularized Least
Squares: Reducing Misspecification Bias with a Flexible and
Interpretable Machine Learning Approach. *Political Analysis*,
22(2):143-168.

## See also

Useful links:

- <https://cetialphafive.github.io/roadrunner/>

- <https://github.com/CetiAlphaFive/roadrunner>

- Report bugs at <https://github.com/CetiAlphaFive/roadrunner/issues>

## Author

**Maintainer**: Jack T. Rametta <jtrametta@gmail.com>
([ORCID](https://orcid.org/0000-0002-9841-146X))

Authors:

- Jack T. Rametta <jtrametta@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9841-146X))
