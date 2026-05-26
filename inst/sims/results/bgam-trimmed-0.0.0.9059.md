# bgam vs baselines — trimmed sim (v0.0.0.9059)

**Date**: 2026-05-26
**Setup**: 4 cells × 50 reps, gaussian DGP1 only (`y = sin(2πx1) + 0.5·x2² + N(0, 0.5²)`). bgam at `mstop=100, nu=0.1, autotune=FALSE`. ares at defaults. ols at full linear model.
**Wall**: 105.5s serial.
**RDS**: `inst/sims/results/bgam-trimmed-0.0.0.9059.rds`

## Test RMSE (mean ± MCSE)

|    n |  p | bgam            | ares            | ols             |
|-----:|---:|----------------:|----------------:|----------------:|
|  200 | 10 | 0.557 ± 0.004   | 0.554 ± 0.004   | 0.854 ± 0.005   |
|  200 | 50 | **0.575** ± 0.005 | 0.617 ± 0.009   | 0.967 ± 0.008   |
| 1000 | 10 | 0.536 ± 0.002   | **0.512** ± 0.002 | 0.841 ± 0.002   |
| 1000 | 50 | 0.538 ± 0.002   | **0.524** ± 0.002 | 0.858 ± 0.002   |

## Mean fit time (seconds)

|    n |  p |  bgam  |  ares  |  ols   |
|-----:|---:|-------:|-------:|-------:|
|  200 | 10 | 0.061  | 0.009  | 0.004  |
|  200 | 50 | 0.277  | 0.087  | 0.007  |
| 1000 | 10 | 0.271  | 0.013  | 0.003  |
| 1000 | 50 | 1.311  | 0.049  | 0.015  |

## Read

- bgam tracks ares within MCSE on small-n, high-p (the niche the spec promised: implicit variable selection on smooth additive signals). Slightly worse than ares on large-n cells where MARS already shines.
- Both crush misspecified `ols` ~30–60% lower RMSE — expected, ols has no smooth basis.
- bgam ~5–25× slower than ares per fit. Acceptable for v1; no autotune in this sim.

## Caveat

Trimmed sim. Full 16-cell × 200-rep run with DGP2 (heteroskedastic) and DGP3 (binomial) deferred to follow-up issue.
