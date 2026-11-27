# Pre-registration: minimal-starter window experiment

Written: 2026-11-27T17:48:48+00:00
predictions.csv SHA-256: `2c7c986d04775ac4f01fbce29fa37889d8553387f0576297ca8fef5fd3c27f5b`
Census used for calibration: `gap_growth_1e8.csv` SHA-256 `f7f8711e7f858c208eb18de2825fd2ad95cda39a6f7f59b5cb0bbe624c3365d5`

Commit this file and predictions.csv BEFORE running run_grid.sh.

## Design
- Gaps: [6, 30, 210, 2310, 30030, 510510]
- Heights (window start): 1e+09, 1e+12, 1e+15, 1e+18, 4e+18
- Window sizes W (even N) x replicates: 1e+05 x 8, 1e+06 x 4, 1e+07 x 2
- Windows: 420; even N examined: 744,000,000
- Engine: engine/goldbach_v6_4_gap.cpp built with -DGAP=g; p_anchor_limit = 1000000; statistic = max minimal g-starter partner over the window (one minpstats row per window).

## Frozen laws: size(M) = C x count, C calibrated on the census (X = 1e4..1e8), never refitted
- size K = exact number of g-starters <= M; size L = S(g) Li2(M)
- count NW = ln N ln W (heuristic); count NN = (ln N)^2 (null: no window-size effect)
- Primary law: K-NW (the heuristic with no approximation in the starter count)

Census residuals (geometric mean of actual/predicted M, by X), constants as frozen:

| law | C | 1e4 | 1e5 | 1e6 | 1e7 | 1e8 |
|---|---|---|---|---|---|---|
| K-NW | 0.599655 | 0.803 | 0.937 | 1.036 | 1.087 | 1.134 |
| L-NW | 0.594783 | 1.084 | 1.007 | 0.993 | 0.973 | 0.957 |
| K-NN | 0.567659 | 0.781 | 0.922 | 1.036 | 1.107 | 1.156 |
| L-NN | 0.563047 | 1.039 | 0.995 | 0.999 | 0.989 | 0.980 |

80% bands come from each law's census residuals.

## Hypotheses and success criteria (decided now, before any window is run)
1. **Window size matters (H1):** each NW law has a smaller mean |log error| than its NN partner (K-NW vs K-NN, L-NW vs L-NN).
2. **Drift continues (H2):** for the primary law K-NW, the geometric-mean actual/predicted rises with height: the value at 4e18 exceeds the value at 1e9.
3. **Accuracy:** report each law's geometric-mean actual/predicted, mean |log error| and band coverage; a law counts as adequate if its overall ratio lies in 0.85-1.15 and 65-95% of windows fall in its band.
4. **Gap uniformity:** for the best law, no gap has a geometric-mean ratio outside 0.80-1.25.

Each hypothesis is reported as SUPPORTED or NOT SUPPORTED as stated, whatever the outcome.
Any later model is reported separately as exploratory, not as a test.
