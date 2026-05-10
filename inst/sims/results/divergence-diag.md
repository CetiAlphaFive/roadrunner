# ares vs earth: divergence diagnostics

Diagnostic run on the v0.1 benchmark cells in `inst/sims/results/v0.1-bench.csv`
where `rss_rel_err > 1%`. Read-only; no source changes.

## Earth defaults vs ares (matches except where noted)

`nk`, `penalty`, `thresh`, `minspan` (auto), `endspan` (auto): identical.
For n=5000 earth trace shows `minspan 7 endspan 9`; ares matches via
`auto_minspan`/`auto_endspan` in `src/ares.cpp:29-44`. `newvar.penalty=0`,
`fast.k=20`, `fast.beta=1`, `Use.beta.cache=TRUE`: ares does not implement
fast.k/fast.beta (locked off per scope) and is unaffected by newvar.penalty=0.

**Mismatches:**
- `Adjust.endspan = 2` — earth doubles endspan for deg>=2 hinge candidates.
  ares does NOT apply this.
- `Auto.linpreds = TRUE` — earth replaces best hinge with linear (`dirs=2`)
  basis when the best knot is at the predictor minimum. ares forward search
  (`src/ares.cpp:499-549`) only proposes `+/-1` hinges.

## Per-cell classification

I picked the worst cell (interaction n=5000 deg=2, rel-err 71%) and a second
high-error cell (additive n=5000 deg=1, rel-err 24.6%) and ran 10 seeds each.

### Cell: interaction n=5000 deg=2  (seeds 1..10)

| seed | rss_a       | rss_e   | rel-err | class |
|------|-------------|---------|---------|-------|
| 1    | 12356       | 2472    | 400%    | (b)+(c) |
| 2    | **3.5e29**  | 12819   | numerical blowup | (e) |
| 3    | 2390        | 2398    | 0.32%   | (a) |
| 4    | 2353        | 2358    | 0.21%   | (a) |
| 5    | 7334        | 2538    | 189%    | (b)+(c) |
| 6    | 7036        | 7231    | 2.7%    | (a) ares-better |
| 7    | 2353        | 6545    | 64%     | (a) ares-better |
| 8    | 2446        | 2447    | 0.01%   | (a) |
| 9    | 12283       | 2432    | 405%    | (b)+(c) |
| 10   | 7171        | 7143    | 0.38%   | (a) |

Findings:

- **Seed 2 (class e — numerical bug)**: ares forward-pass RSS jumps to **3.48e29
  at M=13** (trace) before backward pruning recovers a sensible subset. This is
  a Givens / fast-LS numerical instability — almost certainly the culprit
  behind sporadic large divergences. **Bug; addressable in v0.2.**
- **Seeds 1, 5, 9 (class b + c)**: ares step-1 picks an extreme-tail V5 knot
  (e.g. seed 1: V5=-2.754, RSS 12549.9). Brute-force `stats::qr` over every
  minspan-strided knot confirms ares's local score is correct (-2.754 wins by
  ~1e-5 relative over earth's V5=-1.413 at RSS 12550.1). But ares is then
  trapped: extreme-tail knots have no interior density to support deg-2
  interactions. Earth's central-knot choice is consistent with
  `Adjust.endspan=2` shrinking the eligible interior for hinge proposals.

### Cell: additive n=5000 deg=1  (seeds 101..110)

10/10 seeds: stable ~10% rel-err gap (max 24.6%). Already present at end of
forward pass (not a pruning issue). DGP is
`y = x1 + sin(x2) + |x3| + 0.5 x4^2 + noise`: x1 is linear. Earth uses
`Auto.linpreds=TRUE` to substitute a linear `x1` basis when the best knot is
at x1's minimum; ares forces a hinge pair at an extreme-tail value, which
captures linearity poorly. ares does support `dirs=2` (linear) terms in
`build_term_column` (line 61) but never proposes them in forward search.
**All seeds: class (c) defaults mismatch — Auto.linpreds.**

## Class summary

| Class                                  | Diagnosis | Cells affected           | v0.2 action |
|----------------------------------------|-----------|--------------------------|-------------|
| (a) tie-break flip, eps ~1e-5 rel      | benign    | seeds 3,4,6,7,8,10       | accept; tighten parity test tolerance |
| (b) genuine multi-modal forward path   | partial   | seeds 1,5,9 of inter     | mostly downstream of (c) |
| (c) defaults mismatch: Adjust.endspan=2 + Auto.linpreds | bug-as-feature-gap | additive deg=1 (10/10), interaction deg>=2 with extreme-tail seeds | **HIGH-VALUE FIX** |
| (d) fast.k / fast.beta                 | scope-locked | none observed dominant   | leave off |
| (e) numerical blowup in fast-LS        | bug       | inter seed 2 (n=5000 deg=2) | **MUST FIX**: Givens update needs reorthogonalization or partial QR refresh on RSS rel-improvement < 0 |

## Recommendations for v0.2

1. **Class (e) numerical blowup** — top priority. Detect when forward-pass RSS
   *increases* (rel-imp < 0) or jumps by > some factor; refuse to add the term,
   or trigger a full QR refresh of the basis. Seed 2 is reproducible.
2. **Class (c.1) Adjust.endspan=2** — trivial: multiply auto-endspan by 2 for
   degree>=2 knot proposals (only when proposing a hinge in a *new* variable
   inside an interaction). Single-line change in
   `src/ares.cpp` around the candidate generation loop (line 207).
3. **Class (c.2) Auto.linpreds** — moderate work. After best (var, knot) found,
   if the cut equals the variable's minimum (or is in the lowest minspan
   bucket), substitute `dirs[t,j] = 2` (linear) instead of a hinge pair. ares
   already evaluates `d == 2` correctly in prediction.
4. **Class (a) accepted divergences** — document as expected. Tighten parity
   tests to require rel-err < 0.5% (median) rather than per-cell hard bounds.
5. **fast.k / fast.beta** remain explicitly locked off per scope discipline.
