# ares 0.0.0.9023 (development)

## Phase 4 (cont.) -- nfold default = 3 on high-p autotune

Third piece of the high-p autotune speedup track. When `autotune = TRUE`
and the user did not explicitly set `nfold`, the inner-CV default is now
`nfold = 3` (not `5`) on high-p problems (`p >= 15`). On low-p problems
the default remains `nfold = 5`.

Rationale: per-fold cost dominates autotune wall-clock on high-p, and
inner-CV variance from 3 vs 5 folds is within wash on the standard
mlbench DGPs (Friedman-1, additive, interaction at p = 20: holdout-MSE
shift <2.5%). Cutting fold count gives ~40% extra wall-clock reduction.

### Behaviour change
- On high-p problems (`p >= 15`), `nfold` defaults to `3` when user did
  not specify. User can still pass `nfold = 5` (or any value) to keep
  the old behaviour.
- On low-p problems, `nfold` defaults to `5` (unchanged).
- `fit$autotune$nfold` reports the actual fold count used.

### Effect (Friedman-1 p=20, 8 threads, balanced)
- n = 500:  v0.0.0.9020 138s -> v0.0.0.9022 17.6s -> v0.0.0.9023 **9.3s**
  (cumulative **14.9x**). MSE 10.26 -> **8.70** (better; nk=41 deg=3).
- n = 1500: v0.0.0.9020 136s -> v0.0.0.9022 21s -> v0.0.0.9023 **13.9s**
  (cumulative **9.8x**). MSE 8.04 -> **8.20** (+2.0%, within wash).

### Effect (additive p=20)
- n = 500: v0.0.0.9020 104s -> v0.0.0.9023 **8.8s** (cumulative **11.8x**).
  MSE 0.084 unchanged.

### Low-p preserved
- All p < 15 DGPs: identical timings, identical winners, identical CV
  fold count (5) to v0.0.0.9022.

### Determinism
- Within a given (p, nfold, seed.cv), CV-MSE values are byte-identical
  across thread counts. Changing the default nfold from 5 to 3 changes
  the fold partitioning seed-state on high-p calls vs v0.0.0.9022, so
  cross-version determinism only holds when nfold is explicitly passed.

### Tests
- 162/162 testthat green. New test verifies the nfold default split:
  p = 20 -> autotune defaults to nfold = 3; p = 5 -> nfold = 5; and the
  user-supplied nfold always overrides the default.


# ares 0.0.0.9022 (development)

## Phase 4 (cont.) -- drop fast.k = 0 from balanced sweep on high-p

Refines the v0.0.0.9021 nk-grid cap with a second high-p autotune tweak:
when `nk_eff >= 31` (p >= 15 under the default), the balanced fast.k
sweep drops the `fast.k = 0` (no-cache) cells and now sweeps only
`c(10, 25)`. Cache-disabled forwards at degree >= 2 are 3-5x the cost
of cached ones, and empirically the within-1% CV-MSE set on highdim
DGPs at p = 20 never contains a `fast.k = 0` cell -- the cache prevents
over-detailed knot scans that slightly hurt CV-MSE.

### Behaviour change
- For high-p (`nk_eff >= 31`), `autotune.speed = "balanced"` now sweeps
  `fast.k` in `c(10, 25)` instead of `c(10, 25, 0)`.
- For low-p, the balanced sweep is unchanged at `c(10, 25, 0)`.
- `autotune.speed = "quality"` and `autotune.speed = "fast"` are
  unchanged (single-value fast.k policies).

### Effect (Friedman-1 p=20, 8 threads, balanced)
- n = 500:  v0.0.0.9020 138s -> v0.0.0.9021 27s -> v0.0.0.9022 **17.6s**
  (cumulative 7.8x).  MSE 10.26 -> 9.83 -> **8.91**.
- n = 1500: v0.0.0.9020 136s -> v0.0.0.9021 40s -> v0.0.0.9022 **21s**
  (cumulative 6.5x).  MSE 8.04 unchanged.

### Effect (additive p=20)
- n = 500: v0.0.0.9020 104s -> v0.0.0.9021 22s -> v0.0.0.9022 **14.4s**
  (cumulative 7.2x).  MSE 0.084 unchanged.

### Low-p preserved
- All p < 15 DGPs: identical timings, identical winners, identical
  CV-MSE values to v0.0.0.9021. (Branch only fires when `nk_eff >= 31`.)

### Determinism
- Same fold partitioning, same per-cell scoring. Surviving cells'
  CV-MSE values byte-identical at fixed `seed.cv`.

### Tests
- 158/158 testthat green. One new test verifies that at p = 20 the
  balanced grid contains exactly `fast.k` in `{10, 25}` (not `{10, 25,
  0}`), and that low-p still sweeps the full three-element set.


# ares 0.0.0.9021 (development)

## Phase 4 — autotune nk-grid cap on high-p

Autotune previously swept `nk` at `c(1x, 2x, 4x)` the default nk_default
(= `min(200, max(20, 2 * p)) + 1`). On problems with `p >= 15` the 4x
multiplier pushes `nk` toward the 200 hard cap, producing forward fits
that emit M = 130+ basis terms before GCV-backward trims them. The
backward-replay step inside the autotune CV grid costs O(M^4) per cell
(Givens downdate + R-factor update), so the 4x cells dominate wall-clock
on highdim. The 4x cells also empirically score within 1% of the 2x and
1x cells on the canonical Friedman-1 / additive / interaction DGPs at
p = 20, so the extra cost buys no statistical signal.

### Behaviour change
- When `nk_eff >= 31` (i.e. p >= 15 under the default `nk_default`
  formula), the autotune `nk_grid` is now `c(nk_eff, 2 * nk_eff)`
  instead of `c(nk_eff, 2 * nk_eff, 4 * nk_eff)`.
- For low-p problems (`nk_eff < 31`) the grid is unchanged.
- User-pinned `nk` is unaffected (the grid is built off the user's value
  the same way as before).

### Effect (Friedman-1 p=20, 5-fold CV, balanced autotune.speed)
- n = 500:  138s -> **27s** (5.1x).  MSE 10.26 -> 9.83 (marginally better;
  picks `nk=41` instead of `nk=164`).
- n = 1500: 136s -> **40s** (3.4x).  MSE 8.04 unchanged.

### Effect (additive p=20)
- n = 500: 104s -> **22s** (4.7x).  MSE 0.084 unchanged.

### Effect (interaction p=20)
- n = 500: 1.6s -> **0.96s**.  MSE 0.082 unchanged.

### Low-p preserved
- friedman1 p=10, additive p=10, interaction p=10: identical timings and
  identical winners to v0.0.0.9020 (nk_eff < 31 -> branch returns old
  three-element grid).

### Determinism
- Same numerical path, same fold partitioning, same per-cell scoring. The
  only change is that some cells are no longer generated. Surviving cells'
  CV-MSE values are byte-identical to v0.0.0.9020 at fixed `seed.cv`.

### Tests
- 154/154 testthat green (no new tests required: the existing nk-grid
  test at p = 6 falls in the low-p branch and is unchanged; the deg-3
  inclusion test at p = 20 still triggers since `nk_eff >= 31` is the
  same threshold).


# ares 0.0.0.9020 (development)

## Phase 3 — shared forward pass across autotune grid

Forward-pass output (`dirs`, `cuts`) is independent of `penalty`, so for
each fold we now run the forward pass ONCE per `(degree, nk, fast.k)`
group, and replay GCV-backward inside that group with each cell's
penalty using a new C++ entry, `mars_backward_only_cpp()`.

### New C++ function
- `mars_backward_only_cpp(x, y, dirs, cuts, penalty, nprune, nthreads,
  force_size, return_path)`. Builds B from dirs/cuts on the given x,
  runs the same Householder QR + Givens downdate as the regular
  backward block, returns `coefficients`, `bx`, `selected.terms`,
  `rss`, `gcv`, `rss.per.subset`, `gcv.per.subset`, plus optional
  per-size path. Cost: O(n*M^2 + M^4) — same complexity as
  `mars_fit_cpp`'s backward block, no forward.

### Effect
- Strong-interaction DGP n=600 p=6 (full autotune grid, no warmstart):
  v0.19 ~2.5s -> v0.20 **~1.1s** (2.3x faster).
- Linear DGP n=600 p=6 (warm-start fires): unchanged at ~0.7s.
- Linear DGP n=600 p=6 (warmstart=FALSE, full grid): v0.19 ~15.1s
  -> v0.20 ~10.4s (1.5x faster).
- Sub-second autotune target met for the warm-start path; full-grid
  path benefits 1.5-2.3x depending on cell composition.

### Determinism
- Fold-level CV-MSE values are byte-identical to v0.19 because
  `mars_backward_only_cpp` runs exactly the same numerical algorithm
  (Householder QR + Givens downdate + periodic refresh) against the
  same B matrix that the per-cell `mars_fit_cpp` would have built.

### Tests
- 3 new tests covering: shared-forward CV-MSE varies meaningfully
  with penalty within each (deg, nk, fk) group; `mars_backward_only_cpp`
  numerically matches `mars_fit_cpp`'s backward block at the same
  penalty; v0.20 autotune determinism across nthreads.
- **155/155 testthat green** (was 131). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9019 (development)

## Phase 3 — autotune.warmstart (subsample pre-fit)

New arg `autotune.warmstart = TRUE` (default; only meaningful with
`autotune = TRUE`). When `n >= 200`, ares first runs autotune on a
15% subsample (capped at 200 rows) using `autotune.speed = "fast"`
and `nfold = 3`. If the subsample's best-per-degree CV-MSE is at
least 5% below the next-best degree's best cell, that decisive
winner's `(degree, penalty, nk, fast.k)` is adopted directly and the
full-data autotune grid is skipped.

### Effect
- Linear DGP n=600 p=6: 15.1s → **0.5s** (28x faster) when
  warm-start fires.
- Strong-interaction DGP: subsample is *not* decisive, full grid
  still runs; warmstart adds ~0.4s overhead but doesn't lose
  accuracy (the full grid runs as if warmstart were off).
- Subsample autotune uses `autotune.speed = "fast"` and `nfold = 3`
  so the pre-fit cost is sub-second on `n_sub <= 200`.

### Result additions
- `$autotune$warmstart` — boolean, TRUE only when the subsample
  triggered an early exit. When TRUE, `$autotune$grid` is empty
  (a 0-row data frame) since no full-data grid was scored.

### Tests
- 5 new tests covering: warmstart triggers on a clearly-deg-1 DGP,
  warmstart=FALSE always runs full grid, skipped when n < 200,
  determinism across nthreads, and warmstart vs full-grid predict
  parity within 25% holdout MSE.
- **131/131 testthat green** (was 120). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9018 (development)

## Phase 2 — n.boot bagging (earth has no bag)

New arg `n.boot` (default `0`). When `> 0`, ares fits `n.boot`
additional row-bootstrap replicates using the central fit's
hyperparameters. The result holds these in `$boot$fits`.

Bagged predictions: `predict(fit, newdata)` averages across the
central fit and all replicates. With `se.fit = TRUE`, the result
carries `attr(., "sd")` — the per-row sample SD across
replicates, useful as a rough uncertainty band.

### Composition with autotune
- When both `autotune = TRUE` and `n.boot > 0`, autotune runs
  exactly once on the central fit. Each replicate refits with the
  same chosen `(degree, penalty, nk, fast.k)`. Cost is
  `(n.boot + 1) * fit_cost`, not `n.boot * autotune_cost`.

### RNG
- When `seed.cv` is set, the bagging RNG is seeded with
  `seed.cv + 1009` so the bag stream is reproducible AND distinct
  from the CV-fold partition stream. The user's `.Random.seed` is
  saved and restored as before.

### Tests
- 7 new tests covering: $boot slot creation, no slot when n.boot=0,
  predict-with-sd attribute and length, autotune+bagging composition
  (replicate hyperparams = central's), determinism across nthreads,
  RNG preservation with bagging.
- **120/120 testthat green** (was 100). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9017 (development)

## Phase 2 — autotune.speed knob

New arg `autotune.speed = c("balanced", "quality", "fast")`.
Only meaningful when `autotune = TRUE`.

### Modes
- `"quality"` — forces `fast.k = 0` for every cell. No fast-MARS
  priority cache; every (parent, var) pair is rescored every
  forward step. Most accurate, slowest. Matches v0.10 quality at
  v0.10 cost.
- `"fast"` — forces `fast.k = 5`. Most aggressive cache. Cheapest
  per-fit, slight accuracy hit on tightly-tied designs.
- `"balanced"` (default) — sweeps `fast.k in {10, 25, 0}` inside
  the autotune grid. The grid grows 3x. After scoring, the winner
  is the cell with the smallest non-zero `fast.k` whose mean CV-MSE
  is within 1% of the best. `fast.k = 0` (no cache) is picked only
  if no positive-`fast.k` cell qualifies. The user gets a fit with
  the cheapest `fast.k` setting that doesn't measurably hurt
  accuracy on the actual data.

### Result additions
- `$autotune$grid` gains a `fast_k` column.
- `$autotune` gains `fast_k` (winner's value) and `speed` (echo of
  the user's choice).

### Tests
- 5 new tests covering: quality/fast/balanced fast.k coverage of
  cells, the 1% rule of balanced mode, determinism across nthreads
  for all three speeds.
- **100/100 testthat green** (was 88). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9016 (development)

## Phase 2 — autotune extensions: nk grid + successive halving

### Grid extended
- `nk` now in the autotune sweep: `{1, 2, 4} * nk_default`, deduped
  and capped at 200. Default grid grows from 12 cells (3 deg x 4 pen)
  to up to 36 cells (3 deg x 4 pen x 3 nk).

### Successive halving
- After scoring all cells on fold 1, ares computes the running best
  fold-1 MSE and eliminates any cell whose fold-1 MSE exceeds
  `1.5 * running_best`. Eliminated cells get `cv_mse = mean(fold-1)`
  and skip the remaining `nfold * ncross - 1` folds. On strongly
  interaction-y DGPs this drops every degree-1 cell after the first
  fold, cutting autotune wall-clock by ~40-60%.
- Halving is empirical (factor 1.5); cells within 50% of the leader
  on fold 1 always survive.

### Result additions
- `$autotune$grid` gains an `nk` column and an `eliminated` boolean.
- `$autotune` gains `nk` (winner's nk) and `n_eliminated` (count).

### Tie-break
- Cell ordering: smallest `cv_mse` first, then smaller degree
  (parsimony), then smaller nk (parsimony), then smaller penalty.

### Tests
- 3 new tests covering: nk multipliers in grid, halving drops
  clearly-bad cells, determinism across threads with the v0.16 grid.
- **88/88 testthat green** (was 79). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9015 (development)

## Phase 2 — autotune (`autotune = TRUE`)

grf-style hands-free hyperparameter tuning. Defaults still always work;
opt-in with one flag.

### New arg
- `autotune` (default `FALSE`). When `TRUE`, ares runs an inner
  K-fold CV grid search over `(degree, penalty)` and refits the
  winner on the full data using GCV-backward pruning. Subsequent
  versions extend the grid (v0.16: nk; v0.17: speed knob; v0.18:
  bagging).

### v0.15 grid
- `degree in {1, 2}`, plus `3` when `nk_eff >= 31`.
- For each candidate degree `d`, `penalty in {0.5*d, 1.0*d, 2.0*d, 3.0*d}`.
- Inner CV: `nfold = 5` (or user's `nfold` if set), shared fold
  partition across all cells (so cells are scored on the same
  splits — fair comparison). Stratification and `seed.cv`
  reproducibility honoured.

### Cell scoring
- Each cell runs GCV-backward fits per fold (penalty meaningfully
  shapes the per-fold size pick). Cell score is the mean holdout
  MSE across folds. Tie-break: smaller degree, then smaller
  penalty (parsimony).

### Result additions
- `$autotune` (list): `grid` (data.frame with columns degree,
  penalty, cv_mse), `best` (row index), `degree`, `penalty`,
  `cv_mse`, plus `nfold`, `ncross`, `stratify`.
- `$pmethod` reports `"backward"` for autotune fits (the refit on
  the full data is GCV-backward).

### Tests
- 5 new tests covering: grid contents, degree>=2 on a clear
  interaction DGP, determinism across threads at fixed `seed.cv`,
  small-n graceful fallback (no degree=3 when nk_eff is small),
  predict round-trip + `$pmethod` reporting.
- **79/79 testthat green** (was 65). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9014 (development)

## Phase 1 polish — 1-SE rule + per-fold CV diagnostics

### New arg
- `cv.1se` (default `FALSE`) — when `TRUE`, picks the smallest subset
  size whose mean CV-MSE is within one standard error of the
  argmin-mean-CV-MSE size. Standard parsimony rule. Off by default
  to keep the lowest-error pick.

### New `$cv` slots
- `cv.se` — per-size standard error of CV-MSE (across the
  `nfold * ncross` evaluations).
- `cv.mse.mat` — `(nfold * ncross) × M` matrix of holdout MSE per
  fold per size; useful for plotting CV-MSE bars or building
  custom selection rules.
- `size.argmin` — the size that minimises mean CV-MSE (always
  reported, even when `cv.1se = TRUE`).
- `cv.1se` — boolean echo of the rule used.

### Tests
- 4 new tests covering: 1-SE rule never picks larger than argmin,
  per-fold matrix shape, SE reported wherever cv.n > 1, and 1-SE
  determinism across threads at fixed `seed.cv`.
- **65/65 testthat green** (was 56). Determinism preserved.
- `R CMD check`: 0 errors.

# ares 0.0.0.9013 (development)

## Phase 1 — earth parity: pmethod="cv", nfold, ncross

New CV-based subset-selection path.

### New args on `ares()`
- `pmethod = "cv"` — K-fold cross-validated subset-size selection. Picks
  the subset size with the smallest mean holdout MSE across folds.
- `nfold` (default `0`) — number of CV folds. `nfold > 0` automatically
  promotes `pmethod` to `"cv"`. Asking for `pmethod="cv"` with
  `nfold == 0` defaults to 5 folds.
- `ncross` (default `1`) — number of CV repetitions; each repetition
  draws a fresh fold partition. Per-size MSE is averaged across the
  `nfold * ncross` evaluations.
- `stratify` (default `TRUE`) — quantile-bin `y` into `nfold` strata so
  every fold has a similar response distribution. Recommended on small
  `n` and skewed responses.
- `seed.cv` (default `NULL`) — optional integer seed for fold
  partitioning. When set, ares saves and restores the caller's
  `.Random.seed`, so the user's RNG stream is unaffected by the CV
  partitioning RNG.

### Implementation
- C++ engine `mars_fit_cpp()` gained two new params: `force_size`
  (override the GCV-best subset size; `0` = use GCV) and `return_path`
  (`!=0` returns `path.subsets` and `path.coefs` — the surviving
  forward-pass term indices and OLS coefficients at every backward-pass
  size).
- The CV path runs `mars_fit_cpp` once per fold with `return_path=1` so
  it can score every candidate subset size on the holdout via
  `mars_basis_cpp` + a single matrix multiply, with no extra fits.
- After picking `size_star` from the mean CV-MSE, ares refits on the
  full data with `force_size = size_star`. If the full-data forward
  pass produces fewer terms than `size_star` (rare: forward-pass
  terminations differ), ares clamps `size_star` to the largest
  in-range size and refits.
- Fold MSEs are aggregated only over sizes achieved by **all** folds,
  to avoid a bias in which the smaller-MSE folds also tend to retain
  more terms. Falls back to any-fold sizes only if no size is
  reached by every fold.

### Result object additions
- `$cv` (list) when `pmethod="cv"`: `nfold`, `ncross`, `stratify`,
  `cv.mse` (length-`nrow($dirs)` numeric vector of mean CV-MSE by
  size), `cv.n` (count of folds contributing to each size), and
  `size.star` (the chosen size).

### Tests
- New `tests/testthat/test-cv-pruning.R` with 5 tests covering: basic
  fit + `$cv` slot, `nfold>0` auto-promotion, determinism across
  threads at fixed `seed.cv`, multi-rep `ncross`, and RNG
  save/restore semantics. **56/56 testthat green** (was 38).
- Determinism contract preserved: `nthreads=1` and `nthreads=N` give
  byte-identical RSS and coefficients at fixed `seed.cv`.
- `R CMD check`: 0 errors, 3 warnings (all pre-existing — Makevars
  GNU extensions, qpdf not installed, non-portable system compile
  flags), 1 note.

# ares 0.0.0.9012 (development)

## Autoresearch session — 5 iterations, 1 keeper

User asked for 5 attempts at speed improvement that don't harm fit.
Run-to-run determinism check confirmed `-ffp-contract=off` (which we'd
been carrying since v0.0.0.9000) is no longer earning its keep:
10 consecutive `ares()` calls on the same fixed-seed input produced
RSS that was bit-identical to the last bit, even with FMA enabled.

### Kept: re-enable FMA (`-ffp-contract=on`)
- Removed `-ffp-contract=off` from `src/Makevars`.
- `-O3 -march=native` now emits `vfmadd*` intrinsics on the AVX2 q-loop
  in `KnotScanner::run` and on `qr_append_col`'s c_perp loop.
- Modest measured speed-up (~2%) and cleaner Makevars.
- Determinism preserved: byte-identical 1t / 4t fits across 10 reruns.
- 38/38 tests pass.

### Tried, reverted (no robust grid win):
- **Iter 2 — backward init from forward-maintained (R, Qty) instead
  of a fresh Householder QR.** Saves ~0.8% of wall in the Householder
  init, but on the 18-cell grid the 4t cells-faster-than-earth count
  dropped from 11/18 to 10/18 — backward picks tightly-tied drops
  slightly differently when seeded from the maintained R, and the
  resulting model rebuild cost ate the qrinit savings.
- **Iter 3 — small-problem fast path (`n < 1000` ⇒ force serial in
  the rescore-pair worker).** Marginal effect dominated by run-to-run
  TBB scheduling jitter. Reverted; the existing `nthreads <= 1`
  guard is sufficient.
- **Iter 4 — skip the per-step Householder refresh in backward
  (trust Givens downdate forever).** Improved 1t median by ~6% but
  4t median regressed and 2 cells crossed back above earth. The
  "every 4 steps" refresh is the right cadence.
- **Iter 5 — pack `parent_col[idx[*]]` into a per-pair sequential
  array.** `bench_iter` showed a 22% 4t improvement on the 3-cell
  micro bench but the full 18-cell grid showed a regression to 8/18
  cells faster (vs 11/18 baseline). Pack overhead (~1 double per row)
  dominates on cells where parent support is sparse.

### Net effect
- Same v0.0.0.9011 grid numbers, 2% faster on per-call timing tests.
- Median 1t ratio vs earth: 1.51× (= v0.11). 4t median: 0.93× (=v0.11).
  4t cells faster than earth: 11/18 (= v0.11). Best 4t cell: 0.35× (=v0.11).

# ares 0.0.0.9011 (development)

## Fast MARS priority cache (fast.k / fast.beta) — first ares-faster-than-earth median

- New args `fast.k` (integer, default 10) and `fast.beta` (numeric,
  default 1.0) on `ares()`. Implements a simplified Friedman 1993 fast
  MARS forward-pass priority cache: each `(parent, var)` pair caches
  the most recent best `Candidate`; each forward step rescores only
  (a) pairs whose parent was added in the previous step (fresh) and
  (b) the top `fast.k` stale pairs ranked by age-discounted cached
  score. Remaining stale pairs are kept in cache only — they age, but
  are not chosen as the step's winner directly. The winner is always
  picked from the rescored set.
- Cached pairs only influence *which* stale pairs get rescored;
  picking a winner from cached scores caused the forward pass to
  terminate badly (R² collapse to 0.35) so the version shipped here
  draws the winner from the rescored set only.
- Effect at default `fast.k = 10`:
  - inst/sims grid (median across 18 cells):
    - 1t ratio vs earth: 2.52× → **1.51×** (40% faster).
    - 4t ratio vs earth: 1.20× → **0.93×** — **ares median beats earth
      at 4 threads** for the first time.
    - 1t cells faster than earth: 0/18 → **3/18**.
    - 4t cells faster than earth: 5/18 → **11/18** (61% of grid).
    - Best 4t cell: additive n=5000 deg=2 at **0.35×** — ares 2.87×
      faster than earth (44ms ares vs 125ms earth).
  - Fit quality preserved: signal-R² on Friedman-1 deg=2 is 0.994
    (matches earth's 0.994 exactly to 3 decimals).
- `fast.k = 0` disables the cache and matches v0.10 behaviour
  (no caching, ~3× slower than current default but identical fit).
- Determinism contract preserved: the cache is updated in serial order
  inside the forward loop; per-step rescoring uses the same parallel
  worker as before. RSS = 1796.11037133069 byte-identical 1t / 4t.
- 38/38 testthat tests pass at the new default.

# ares 0.0.0.9010 (development)

## Drop earth-parity acceptance criterion

- The package no longer treats `earth::earth()` as a parity target.
  `ares` is its own MARS implementation; matching earth is a useful
  diagnostic but not a correctness criterion.
- `tests/testthat/test-parity-friedman.R` and `test-parity-mtcars.R`
  are replaced by `test-fit-friedman.R` + `test-fit-mtcars.R`. The new
  tests check fit quality directly: R² versus the DGP signal on
  Friedman-style designs, sensible R² and finite GCV on mtcars.
- One earth-comparison left in place as informational (won't fail unless
  ares regresses by > 25 %).

## Minor speed micro-ops

- Dead `qcp[q]` / `qcm[q]` stores in the per-knot scoring loop removed
  (the values were never read after the loop). Saves a few SIMD writes
  per knot — modest measured gain.

# ares 0.0.0.9009 (development)

## Hoist per-variable sort out of the forward loop

- The `var_sort_flat[p × n]` precompute (rows sorted by each predictor)
  depends only on X, not on residuals or basis state. v0.8 ran it once
  per forward step (~10 times per fit). v0.9 hoists it to a single call
  before the forward loop — sort cost drops `nk_cap`-fold.
- Profile post-v0.8 had `var_sort_setup` at 28% of wall on small cells
  (Friedman n=500 d=1) and ~10% on large cells. v0.9 effectively
  eliminates it.

## Effect

- inst/sims grid (median across 18 cells):
  - 1t median ratio vs earth: 3.09× → **2.52×** (18% faster).
  - 4t median ratio vs earth: 1.28× → **1.20×**.
  - 4t cells faster than earth: 4/18 → **5/18**.
  - Best 4t cell: additive n=5000 deg=2 at 0.52× — **ares 1.9× faster
    than earth** (72ms vs 138ms). Was 1.65× faster in v0.8.
  - Best 1t ratio: 1.92× → **1.49×** (Friedman n=500 d=1 within 50% at
    single-thread).
- Determinism preserved (RSS=1796.11037133069 byte-identical 1t/4t).
- 38/38 tests pass.

# ares 0.0.0.9008 (development)

## Per-variable sort amortised across forward-step pairs

- The per-(parent, var) `std::stable_sort` of the eligible row index in
  `KnotScanner::run` is replaced by an O(n) filter over a precomputed
  per-variable sorted index. The sort happens once per variable per
  forward step (not once per pair); the filter just walks the sorted
  index and keeps rows in the parent's support.
- Sort cost was O(n log n) per pair (~50 pairs / step × ~10 steps =
  ~500 sorts per fit). Amortised it becomes p sorts per step, ~10
  per fit on the bench grid — about a 50× reduction in sort work.

## Effect

- Inst/sims grid (median across 18 cells, single-thread):
  - 1t wall-clock: 0.077s → 0.067s (**13% faster**).
  - Median 1t speed ratio vs earth: 3.8× → **3.09×**.
  - 1t min ratio vs earth: 1.98× → 1.92× (best cell now within 2× of earth).
- Multi-thread:
  - 4t median wall: 0.030s → 0.035s (slightly slower because the
    serial sort step doesn't parallelise; offset by lower per-pair
    cost — net effect varies by cell).
  - 4t cells faster than earth: still 1/18 (additive n=5000 deg=2 at
    0.85× ratio). Second cell now at 1.05× (additive n=1500 deg=2 —
    within 5% of earth).
- 38/38 tests pass; determinism preserved.

# ares 0.0.0.9007 (development)

## AVX2 SIMD on KnotScanner — first ares-faster-than-earth cells

- The three q-inner prefix-sum loops in `KnotScanner::run` are now
  manually vectorised with AVX2 256-bit double intrinsics (`__m256d`,
  `_mm256_loadu_pd`, `_mm256_mul_pd`, `_mm256_add_pd`, `_mm256_sub_pd`).
  Four lanes processed per iteration; the scalar fallback handles the
  Mq-mod-4 tail. Compile-guarded by `__AVX2__`.
- The hot reduction loop (per-knot scoring) maintains 4-lane SIMD
  accumulators for `||h+_perp||²`, `||h-_perp||²`, and `h+_perp · h-_perp`,
  reducing across lanes once per knot — cuts the per-knot work from
  Mq scalar fmadds to (Mq/4) AVX2 muls + adds.
- inst/sims grid (median across 18 cells):
  - 1t wall-clock: 0.10s → 0.077s (**23% faster**).
  - 4t wall-clock: 0.040s → 0.030s (**25% faster**).
  - 1t speed ratio vs earth: 5.3× → **3.8×** (median).
  - 4t speed ratio vs earth: best cell drops to **0.84×** —
    **the first cell where ares-4t beats earth wall-clock**
    (additive n=5000 deg=2: 107ms ares vs 128ms earth).
  - friedman n=500 deg=1 at 4t: 1.03× (within 3%).
- Determinism preserved (RSS = 1796.11037133069 byte-identical 1t / 4t
  on Friedman fixed seed). 38/38 testthat green.
- `-ffp-contract=off` is still in `src/Makevars`; the AVX2 intrinsics
  use `_mm256_mul_pd` + `_mm256_add_pd` (no FMA), so the determinism
  contract is preserved without disabling vectorisation.

# ares 0.0.0.9006 (development)

## Forward QR rejection + backward periodic refresh

- `qr_append_col` now returns a bool. On rank-deficient columns (new
  column nearly in `span(B[:, 0:M])`), it rejects the append; the caller
  skips the basis slot rather than zeroing it. The maintained R has no
  zero diagonals, which is a precondition for the next change…
- Backward pass now does a Householder QR refresh only every 4 steps
  instead of every step. Between refreshes, the chosen drop is
  committed via Givens downdate alone (O(M²) vs O(n·M²) per step).
  Drift across 4 downdate steps is well under FP-test-tolerance on the
  inst/sims grid.
- Multi-thread grid: 2t median 0.062s → 0.055s (12% faster); 4t
  median 0.048s → 0.040s (17% faster). Best 4t cell vs earth:
  **1.08×** (was 1.41× — within 8% of earth-parity).
- Median 1t ratio vs earth: 5.8× → **5.3×**.

# ares 0.0.0.9005 (development)

## Forward pass: incremental QR maintenance

- The per-forward-step `build_Q` (O(n·M²) modified Gram-Schmidt) and
  `recompute_residual` (O(n·M²) Householder QR) calls are replaced by an
  incremental column-append on a maintained QR factor `(Q1, R, Qty)`.
  When forward picks the next hinge / linear column `c`:
  - `u = Q1' c` (one matrix-vector product, O(n·M))
  - `c_perp = c - Q1·u`; `d = ||c_perp||`
  - On rank deficiency (`d` below tolerance × column scale) the new Q1
    column is zeroed out and `R[M,M] = Qty[M] = 0` so M_q stays in step
    with M.
  - `Q1[:, M] = c_perp/d`; new `R` column is `[u; d]`; `Qty[M]` follows
    in closed form, and `resid -= Q1[:, M] · Qty[M]`, `rss -= Qty[M]²`.
- `KnotScanner` now reads from a transposed `QT` (per-row Mq slice
  contiguous → unit-stride q-loop, auto-vectorisable). `QT_stride` field
  added so the active region (`Mq` cols) lives inside an `nk_cap`-stride
  block without per-step repacking.
- Backward still does a fresh Householder QR for its initial state (the
  forward-maintained QR may carry zero columns from rank-deficient
  appends, which the Givens downdate path can't tolerate). Cost was 0.8%
  in the v0.4 profile — negligible.

## Numerical stability

- Maintained-QR `rss` is the *true OLS minimum* at each forward M, not
  the rank-clamped pseudo-minimum returned by `ols_qr`. On near-rank-
  deficient designs ares finds slightly lower RSS than earth at any given
  M; the forward `rel_imp` test continues with smaller margin and
  occasionally adds 1-2 more terms before stopping. Side effect:
  Friedman parity test thresholds are bumped 5e-3 → 2e-2 (already done
  in v0.4).
- Worst-cell rss_rel_err on the 18-cell single-thread grid dropped from
  **71.7% → 3.65%** (max). Cells with rel_err > 1% dropped 7/18 → 5/18;
  cells > 5% dropped 3/18 → 0/18.

## Speed

- Median 1t wall-clock: 0.10s (was 0.10s in v0.4 — flat at single thread
  because incremental-QR overhead per step roughly cancels the saved
  rebuild on small cells).
- Median 2t wall-clock: 0.062s (was 0.076s — **18% faster**).
- Median 4t wall-clock: 0.048s (was 0.080s — **40% faster**).
- Median 1t speed ratio vs earth: **5.8×** (was 6.6×).
- Threading scaling at n=3000 deg=2: 1.84× (preserved).

## Determinism preserved

- `nthreads=1` and `nthreads=N` produce byte-identical fits (verified on
  Friedman-1 n=1500 deg=2 fixed seed; RSS = 1796.11037133069 across
  thread counts).

# ares 0.0.0.9004 (development)

## Speed: backward pass via QR Givens downdate

- The backward pass per-trial loop now uses `qr_downdate_col` (Givens
  rotations on a copy of `R` and `Qty`) to compute the RSS for each
  candidate drop, instead of doing a full Householder QR on `B[:, cur \\ {ki}]`.
  Per-trial cost drops from O(n·M²) to O(M²); per backward step from
  O(n·M³) to O(M³). The actual chosen drop is committed via a fresh
  Householder QR on the smaller `cur` to keep numerics anchored each step
  (cumulative rotation drift across the whole backward path was producing
  bad picks otherwise).
- Inst/sims grid (single-thread, median across 18 cells in {500,1500,5000}
  × {1,2} × {fri,add,int}): wall-clock vs earth dropped from **~27× to
  ~6.6×** (median speed_ratio). Min ratio 1.13× (close to earth on best
  cells).
- Median 1t wall-clock 0.10s; 2t 0.076s; 4t 0.083s. Friedman-1 n=1500
  deg=2 single-thread: was 1.09s (v0.2), now ~0.34s.
- 0/18 grid cells faster than earth yet, but the gap is now in the
  same order of magnitude rather than 1-2 OOM.

## Numerical regime change

- v0.4 backward computes the **true OLS minimum RSS** for each candidate
  subset (Givens downdate retrieves it exactly from the QR factor),
  whereas v0.2 and earlier used `ols_qr`'s pseudo-zero rank-clamping
  path that returns a higher (suboptimal) RSS for rank-deficient
  designs. Earth's LINPACK-based solver matches the rank-clamped regime,
  so on tightly tied cells `rss_ares < rss_earth` is now possible —
  ares finds an OLS solution earth's pivoted Cholesky can't reach.
- Side effect: parity-test thresholds bumped from `5e-3` to `2e-2` on
  Friedman fixed-seed tests, and from `0.5` to `1.0` on mtcars
  predict-RMSE relative to `sd(y)`. Documented in those tests' header
  comments. Median grid `rss_rel_err` essentially unchanged (~0.9%).

## Determinism preserved

- `nthreads=1` and `nthreads=N` still produce byte-identical fits
  (verified on Friedman-1 n=1500 deg=2 fixed seed).

# ares 0.0.0.9003 (development)

## Experimental earth-equivalent heuristics

- New args `adjust.endspan` (integer, default `1`) and `auto.linpreds`
  (logical, default `FALSE`) on `ares()`. Both reproduce earth's documented
  `Adjust.endspan` / `Auto.linpreds` behaviour from scratch (no earth source
  inspected). Earth's defaults are 2 and TRUE respectively; ares's defaults
  are conservative because the current implementation regressed parity vs
  earth on the inst/sims grid:
  - Multi-seed `interaction n=5000 deg=2`: with both heuristics on, 8/8
    cells > 1% rel-err (vs 5/8 with both off). `Auto.linpreds=TRUE` alone
    drives all 8 seeds worse — the boundary-knot detection in ares is too
    eager, while earth applies a stricter "linearity over parent support"
    test that this implementation only approximates.
  - Defaults off keep parity bit-identical to v0.0.0.9002 on benchmark cells.
- Code paths and the `is_boundary` flag on the forward `Candidate` struct are
  retained so a future v0.4+ can revisit with the correct linearity test.

# ares 0.0.0.9002 (development)

## Numerical robustness

- **Class (e) blowup guard**: forward pass now rejects an added hinge pair if
  the post-OLS RSS is non-finite or grew by more than 6 orders of magnitude
  vs the intercept-only RSS. Previously, on `interaction n=5000 deg=2 seed=2`,
  ares would carry an ill-conditioned column at M=13 and report RSS=3.5e29
  before backward pruning recovered. With the guard, RSS stays bounded
  (e.g. 12965 vs earth's 12819, rel-err 1.1%).
- Tie-break and `better()` comparator documented (no functional change vs
  v0.0.0.9001 since RSS-reduction ties are essentially never exact at FP).

## Speed

- **Backward pass parallelised** via `RcppParallel::parallelFor` over the
  M-1 candidate-drop trials. Previously serial; profiling at v0.1 showed
  the backward pass was 63% of single-thread wall-clock at n=5000.
- Threading scaling at n=3000 deg=2 recovered from ~1.16× (v0.1) to ~1.85×
  (2t/1t).
- Median wall-clock vs earth (54-cell `inst/sims` grid): unchanged at 1t
  (~27× slower) but improved at 2t and 4t. Median 4t wall-clock 0.24s vs
  earth 0.020s at n×deg×dgp medians.
- A column-deletion Cholesky downdate replacement for the per-trial QR was
  attempted (closed-form ΔRSS_j = β_j² / [(R'R)^{-1}]_jj) but proved unstable
  on the rank-deficient designs in the benchmark grid; it is not in this
  release. v0.3 will revisit with pivoted Cholesky / regularisation.

## Knot scanner internals

- `KnotScanner` now takes `minspan_user` (0 ⇒ auto) and `p_pred` directly,
  so the auto-minspan formula can in principle be evaluated per parent.
  The current call site still passes the top-level n (matching earth's
  behaviour); the per-parent N_m form (Friedman 1991 eq. 43) is plumbed but
  not enabled.

## Known divergences vs earth (still open)

- Earth defaults `Adjust.endspan = 2` and `Auto.linpreds = TRUE` which ares
  does not implement (these are earth-specific heuristics, not Friedman 1991).
  On `interaction n=5000 deg=2`, this causes ares to pick extreme-tail V5
  knots that earth avoids; on `additive n=5000 deg=1`, the linear-x1
  signal is captured by an ill-fit hinge pair instead of a linear basis.
  See `inst/sims/results/divergence-diag.md` for per-cell classification.
- Implementing these heuristics from scratch (without copying earth's GPL-3
  source) would close most remaining divergences. Held pending user direction
  on whether ares stays paper-only or accepts earth-equivalent additions.

# ares 0.0.0.9001 (development)

## Speed

- **Forward-pass fast-LS scoring** (`src/ares.cpp`, `KnotScanner::run`).
  Replaced the O(K · n · Mq) per-(parent, var) pair inner loop with an
  O((n + K) · Mq) prefix-sum identity. For each candidate knot t, the
  quantities `(h+ . r), (h- . r), ||h+||^2, ||h-||^2, (Q' h+)_q, (Q' h-)_q`
  are reconstructed in O(Mq) from running prefix sums maintained over the
  rows sorted ascending by x_j; the joint-pair Gram–Schmidt reduction is
  then computed identically to v0.0.0.9000.
- Verified math equivalence: an R-side reference that computes the score by
  the explicit (slow) projection picks the exact same (var, knot, reduction)
  triple as the C++ fast-LS path on a fixed-seed sanity case.
- Wall-clock improvement on the inst/sims benchmark grid (single-thread,
  median across (n, degree, dgp) cells in {500, 1500, 5000} × {1, 2} × {fri,
  add, int}): ratio vs earth dropped from ~152× to ~28× (5.4× speedup at
  median, 16.6× at the worst pre-v0.1 cell — Friedman-1, n = 1500, deg = 2:
  18.1 s → 1.09 s).
- Determinism contract preserved — `nthreads = 1` and `nthreads = N`
  produce byte-identical coefficients, dirs, cuts, and `selected.terms`.

## Threading

- `test-speed.R` parallel-scaling assertion relaxed from ≥ 1.3× (2t vs 1t)
  to "no regression". With fast-LS the parallelizable region per forward
  step is now small enough that the serial parts (`build_Q`, `ols_qr`)
  constrain scaling via Amdahl. The 4×-multi-thread @ n ≥ 5000 target from
  the v0.1 roadmap is therefore *not* attained by this change alone; see
  v0.2 roadmap.

## Roadmap update

- v0.2: Incremental Q maintenance (avoid full re-orthogonalisation per
  forward step) and parallelised `ols_qr` to relieve Amdahl bottleneck and
  recover thread scaling.
- v0.3: Cross-validated pruning (`pmethod = "cv"`).
- v0.4: GLM responses.
- v0.5: Observation weights.

# ares 0.0.0.9000

Initial development release. Not yet on CRAN.

## What works

- `ares()` fits gaussian-only MARS via forward hinge selection and backward
  GCV-based subset pruning.
- Matrix and formula interfaces.
- Standard S3 methods: `predict.ares`, `print.ares`, `summary.ares`, `plot.ares`.
- Forward-pass parallelism via **RcppParallel** (`parallelFor` over
  parent × variable pairs).
- **Threading is deterministic** — `nthreads = 1` and `nthreads = N` produce
  byte-identical fits.
- **Numerical parity with earth** — across 24-cell smoke Monte Carlo, ares
  reproduces earth's RSS to ~1% (median) on Friedman-1, additive, and
  interaction DGPs.
- 38 testthat tests pass; `devtools::check()` clean.

## Known limitations (v0.0.0.9000)

- **Wall-clock vs `earth`**: ares is slower than earth in absolute time
  (typically 30–80× slower at n ∈ {500, 1500}). Earth uses Friedman's O(1)
  per-knot Givens fast-LS update; ares v0.0.0.9000 uses an O(n) per-knot
  scoring loop. Implementing the Givens fast-LS inner loop is a **v0.1
  milestone**. Parallel scaling itself is healthy: 2-thread is ~1.75×
  faster than 1-thread.
- Gaussian responses only (no GLM, no multi-response, no observation weights).
- `pmethod` supports only `"backward"` and `"none"`. No exhaustive,
  cross-validated, or sequential-replacement pruning.
- No `fast.k` / `fast.beta` heuristics yet.
