# Pre-registration, Round 2: worst-N residue correction for extreme starter witnesses

Written: 2026-11-27T18:19:48+00:00
predictions_r2.csv SHA-256: `f4df2bff5d8d2197fbdb74ba2a7f4ab2e9f849017499acddc81370d423037b3b`
Calibration input: the 420 Round-1 windows, minpstats files in `runs/` (combined SHA-256 `a05ee4a108f9b5178a03f72cd687d5ead19cc24f1e97e626942890b40fe6752d`),
with window metadata from `prereg/predictions.csv` (SHA-256 `2c7c986d04775ac4f01fbce29fa37889d8553387f0576297ca8fef5fd3c27f5b`).

Commit this file and predictions_r2.csv, and push, BEFORE running run_grid_r2.sh.

## Design
- Gaps: [6, 12, 42, 66, 330, 390, 2310, 4620] (6 and 2310 carried over from Round 1; 12 and 4620 share their odd prime divisors)
- Heights (window start): 3e+09, 3e+12, 3e+15, 2e+18, 8e+18 (none used in Round 1)
- Window sizes W x replicates: 1e+05 x 8, 1e+06 x 4, 1e+07 x 2
- Windows: 560; even N examined: 992,000,000
- Engine: engine/goldbach_v6_4_gap.cpp with -DGAP=g, p_anchor_limit = 1000000; results in runs_r2/.

## Frozen laws (fitted on Round 1, never refitted)
F(g) = product over primes 5 <= l <= 10^5 with l not dividing g of (1 - 1/(l-2)^2).

| law | C | ln N exponent | Round-1 actual/pred | Round-1 mean abs log error |
|---|---|---|---|---|
| R2-F | 0.645853 | 1.0000 | 0.996 | 0.109 |
| R2-FH | 0.431996 | 1.1157 | 0.995 | 0.103 |
| K-NW | 0.688266 | 1.0000 | 0.999 | 0.118 |

Primary law: R2-F. 80% bands come from each law's Round-1 residuals.

## Mechanism statistics (frozen)
- Unfavourable class: N != 0 (mod l), and N != -g (mod l) when l does not divide g. Primes tested: [3, 5, 7, 11, 13].
- H_g(N) = product over odd primes l <= 97 of (1 - b/l)/((1 - a/l)(1 - 1/l)),
  a = |{0, -g} mod l|, b = |{0, -g, N} mod l|.
- For each window: the hardest N (minpstats N_at_max) is ranked against 1000 random even N from the same
  window, drawn with numpy default_rng seeded by the first 16 hex digits of SHA-256("r2:g:start:end").
  Mid-rank = (#sample with H < H_hard + 0.5 * #ties) / 1000.

## Hypotheses and criteria (decided now; each reported as SUPPORTED / NOT SUPPORTED)
- **H1 (gap correction helps):** R2-F has a smaller mean |log error| than K-NW over all Round-2 windows.
- **H2 (height law holds):** for R2-FH, the geometric-mean actual/predicted at 8e18 divided by that at 3e9 lies in 0.90-1.10.
- **H3 (same divisors, same behaviour):** under R2-F, the geometric-mean actual/predicted of g=12 over g=6, and of g=4620 over g=2310, each lie in 0.90-1.10.
- **H4 (hardest N are unfavourable):** for each l in [3, 5, 7, 11, 13], the number of windows whose hardest N is unfavourable at l
  exceeds its expectation under random even N, one-sided normal approximation to the Poisson-binomial, p < 0.001.
  Supported if it holds for every l.
- **H5 (hardest N have small H):** the mean mid-rank of H_g(hardest N) is below 0.40, and more than half the windows have
  mid-rank below 0.5 (one-sided sign test, p < 0.001).
- **Accuracy (reported for each law):** adequate if overall actual/predicted lies in 0.85-1.15 and 65-95% of windows fall in its band.

No further laws are added after registration; anything else is reported as exploratory.
