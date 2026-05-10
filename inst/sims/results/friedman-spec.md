# Friedman MARS spec vs ares — discrepancy analysis

**Sources** (open mirrors, 2026-05-10):

- Friedman 1991, "MARS," *Annals of Statistics* 19(1):1–67 — `stat.yale.edu/~lc436/08Spring665/Mars_Friedman_91.pdf`. Equations cited by the paper's numbering.
- Friedman 1993, "Fast MARS," Stanford TR-110. Public mirror PDF is image-only (no OCR available in this environment). Supplementary detail from milbo's `earth-notes.pdf`.
- ares: `/home/jack/Dropbox/ares/src/ares.cpp` (800 lines).

## Forward-pass scoring

**Paper.** Algorithm 2 (§3.4, p.17). For each parent `B_m`, variable `v`, knot `t ∈ {x_{v,j} | B_m(x_j) > 0}`:

```
g ← Σ a_i B_i(x) + a_M B_m(x)[+(x_v − t)]_+ + a_{M+1} B_m(x)[−(x_v − t)]_+
lof ← min_{a_1..a_{M+1}} LOF(g) = Σ (y_i − g(x_i))^2     (eq. 29)
```

The paper jointly refits **every** coefficient `a_1..a_{M+1}` at every candidate; the Cholesky-update machinery of §3.9 makes that affordable.

**ares.** `KnotScanner::run` (lines 247–391) does not full-refit. It scores the joint RSS reduction by projecting the current OLS residual `r = y − Bβ̂` onto `span(h+·B_m, h−·B_m)` after Gram–Schmidt orthogonalization against the precomputed orthonormal basis `Q` of `B`:

```
red_p = (cp·r)^2 / ||cp_⊥||^2
red_m = (e2·r)^2 / ||e2||^2     (2nd GS step against cp_⊥)
RSS_red = red_p + red_m
```

Algebraically equivalent to the LOF reduction in eq. 29 because `r ⊥ span(B)` after `recompute_residual`. Differs only in numerical path.

## Knot enumeration order

**Paper.** §3.9 (p.29) explicitly says the inner loop visits knots **in decreasing order**: *"Friedman and Silverman (1989) developed updating formulae for least-squares fitting of q = 1 splines by visiting the eligible knot locations in decreasing order…"* The eq. 52 rank-1 update assumes each new `t` is the next-smaller value.

**ares.** `KnotScanner` sorts *ascending* in `x_j` (line 256) and sweeps `k = 0…n_eli-1`. Prefix sums `S_*` are low-side; high-side reconstructed as `H_* = T_* − S_* − row_k`. Mathematically equivalent — prefix-sum identities are direction-agnostic — but FP summation order differs.

**Practical consequence:** in exact arithmetic, none. Under FP, ULP-scale differences in `red_p + red_m` can flip near-tied candidates. Matters mainly via tie-break (next section).

## Tie-break

**Paper.** No explicit tie-break. Algorithm 2 line 8 uses *strict* `if lof < lof*`, so the **first candidate** in the loop order wins. With `m` outer, `v` middle, `t` inner *descending*, "first" = largest `t`, smallest `m`, smallest `v`.

**ares.** `better()` + inline (lines 207–213, 362–369): largest `rss_red`; tie → smallest `var`; tie → **smallest** `cut`; tie → smallest parent.

**Discrepancy.** ares prefers smallest `t` on a tie; the paper's loop order (descending) implies largest `t`. Likely the dominant cause of forward-pass divergence from earth on near-tied knots.

## GCV / pruning

**Paper.** §3.6 eq. 30–32:

```
GCV(M) = (1/N) Σ (y_i − f̂_M(x_i))^2 / (1 − C̃(M)/N)^2
C̃(M) = C(M) + d · M             (eq. 32)
C(M) = trace(B(B^T B)^{-1} B^T) + 1 = M + 1     (eq. 31, the linearly-indep count)
```

Default `d = 3` for general MARS, automatically reduced to `2d/3 ≈ 2` for additive (degree=1) models (§4, p.32).

**ares.** `compute_gcv` lambda (lines 651–656):
```
C = M + penalty * (M − 1) / 2
GCV = RSS / (n · (1 − C/n)^2)
```
This is the earth / ESL form (each hinge pair contributes one knot ⇒ `(M-1)/2` knots), not Friedman's literal eq. 32. For parity vs earth this is the **correct** choice. Backward pass (lines 658–706) matches Algorithm 3 verbatim. `2d/3`-on-additive auto-reduction must live in the R wrapper, not C++; verify that default matches earth's `penalty = if (degree > 1) 3 else 2`.

## Fast-LS Cholesky

**Paper.** §3.9 (eq. 47–53). Friedman uses **rank-1 Cholesky update** of the Gram matrix `V = B^T B` and RHS `c = B^T y`. Eq. 52 gives an explicit closed form for the (M+1)-st row/column of `V` at each new `t` (in decreasing-`t` sweep), in `O(M·N_m)` per knot; Cholesky update of the factor is `O(M^2)`. Total inner cost: `O(M·N_m + M^2·N_m/L)` per (parent, var). **Not Givens rotations** — the CLAUDE.md handoff's "Givens fast-LS" wording does not appear in either Friedman paper.

**ares.** Uses **prefix-sum identities over sorted-x** (lines 222–390), not Cholesky update. Every quantity entering `red_p + red_m` is linear in `t` with coefficients that are sums over low-side / high-side rows; sums maintained as running prefix updates at `O(M_q)` per row. Per-knot cost `O(M_q)`; per-pair `O((n+K)·M_q)`.

**Equivalence?** Algebraically yes, given `r ⊥ Q`. Reproduces the same RSS reduction Friedman gets via Cholesky, up to FP summation order. Diverges when (a) `Q` is rank-deficient (ares drops; earth keeps via pivoted Cholesky w/ ridge — eq. 54), (b) cancellation on near-collinear `B_m·x_v` columns, or (c) the `tiny = 1e-12` Gram-Schmidt threshold (line 305) skips a candidate earth would keep.

## Delta list — ranked by likely impact on earth-parity

| # | What | Paper says | ares does | Likely parity impact |
|---|------|------------|-----------|----------------------|
| 1 | Knot tie-break (within a (parent, var)) | implicit: first-encountered in *descending* sweep ⇒ **largest** `t` | smallest `t` | **High.** Different choice on near-ties means different forward-pass tree → different terms → different RSS / GCV at every subsequent step. Most likely root cause of ares-vs-earth term divergence. |
| 2 | Knot sweep direction | descending in `x_j` (§3.9 line 1405) | ascending in `x_j` (line 256) | **Medium.** Equivalent in exact arithmetic; differs by ULP-scale rounding under FP. Matters only via tie-break above. |
| 3 | LOF computation method | rank-1 Cholesky update of `B^T B` (eq. 52) | prefix-sum identities + Gram–Schmidt projection against orthonormal `Q` | **Low–medium.** Algebraically equivalent. ares re-derives Q each step (`build_Q`); earth maintains the Cholesky factor incrementally. Numerical paths differ → ULP-scale RSS diffs. |
| 4 | Rank-deficiency handling | pivoted Cholesky (eq. 54: ridge-like `(V + εD)a = c`) | drop-column with `tol = 1e3 · eps · scale` in `build_Q` (lines 466–469), pseudo-zero in `ols_qr` (line 178) | **Low** unless data is very degenerate; can flip a candidate near rank boundary. |
| 5 | Penalty algebra | `C̃(M) = (M+1) + d·M` with `d = 3` (additive: `d ≈ 2`) | `C = M + p·(M−1)/2`, no auto-additive reduction | **None vs earth.** ares matches *earth's* form, not Friedman's literal eq. 32. Earth is our parity target → correct choice. |
| 6 | Hinge sign order added to basis | paper unspecified; Algorithm 2 lines 9–10 add `[+(·)]_+` then `[−(·)]_+` | adds `s = +1` first, then `s = −1` (line 608) | **None.** Just term-index permutation; final GCV invariant. earth is documented as adding `−` first; if ares strictly tracks earth term-indices for diagnostics there is a small naming-only difference. |
| 7 | Backward pass | Algorithm 3: drop term minimizing LOF among non-intercept terms, refit, track best GCV | identical (lines 658–706) | **None.** |
| 8 | minspan (eq. 43) | `L(α) = ⌊−log₂(−(1/(nN_m))·log(1−α))/2.5⌋` with `α=0.05` | `auto_minspan(p, n) = ⌊−log₂(−(1/(p·n))·log(1−0.05))/2.5⌋` (lines 86–93) | **Low.** Uses fixed `n` and `p` instead of `N_m` (current parent's nonzero count). Friedman writes `nN_m`; ares uses `pn`. For full-support parents (intercept) `N_m = N` and `n` (Friedman's predictor count) = `p` (ares), so they agree. For nested terms `N_m < N`, so ares's minspan is too *small* relative to the paper. earth tracks `N_m` per-parent. |
| 9 | endspan (eq. 45) | `Le(α) = ⌊3 − log₂(α/n)⌋` with `α=0.05` | `auto_endspan(p) = ⌊3 − log₂(0.05/p)⌋` (lines 95–101) | **None.** Matches paper exactly with `α=0.05` and `n → p`. |
| 10 | Hinge basis function | `b(x|s,t) = [s(x−t)]_+` (eq. 33) | `hinge(s,x,t) = max(0, s·(x−t))` (lines 104–107); `s ∈ {+1, −1}` with internal `dirs=2` reserved for linear (unused) | **None.** |
| 11 | Continuity (piecewise cubic) | §3.7: replace each linear hinge with a piecewise cubic of central knot `t` and side knots `t_±` after the model is selected | not implemented (documented in `ARCHITECTURE.md` line 164 as out-of-scope) | **None vs earth's default.** earth defaults to piecewise linear (`Use.beta.cache`-related continuity options off). For users who set earth's `Use.cubic = TRUE` ares will diverge by design. |
| 12 | Fast MARS heuristic (1993 TR) | priority queue / "h" parameter limiting parent re-evaluation | not implemented (CLAUDE.md scope discipline) | **None vs earth's default `fast.k = 20`.** earth defaults *do* use Fast MARS; if v0.1 wants strict parity with earth defaults this becomes relevant. To match earth bit-exact, also disable `fast.k` on the earth side (`fast.k = 0`). |

## Bottom-line ordering for the v0.1 parity push

If "match earth's term selection at FP-tied knots" is the goal, the highest-leverage fix is **#1 (tie-break to prefer largest `t`)** combined with **#2 (descending sweep)** so the comparator order in `KnotScanner` matches earth's conventional first-found-wins behavior. Items #3 and #4 are second-order (ULP-only) and only matter when many candidates are within `~1e-10·RSS` of each other.

Once #1–#2 are aligned, remaining parity gap should be dominated by #4 (rank tolerance) for high-correlation designs and #8 (`minspan` arg) for deeply-nested terms.
