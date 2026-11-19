# Pre-registration: extreme-value test of the Goldbach residual null
# Author: Santhiagu Antony Raj Parthibanathan. Drafted 2026-11-19 (IST), BEFORE the campaign runs.
# Commit this file to the repo (or a Zenodo deposit) before launching the runs below.

## Null hypothesis H (unchanged from the N- and B-scaling tests)
After exact removal of congruence-explained failures (candidate p where N - p has a
prime factor <= B), residual failures of minimal-partner trials at height N are
independent Bernoulli events at the local prime density among B-rough odd numbers.

## Campaign (five runs, v6.4, OCI E5, 48 threads)
Each run: span 10^14 of N, minp-per-window=100 => 10,000 windows x 100 rows = 1,000,000
minimal partners per scale. Settings: block_bits 24, p_anchor_limit 10^7, sample_limit
100000, qhot_size 512.

  R1: [1*10^15,        1.0001*10^15      ]   R2: [1*10^16, 1.0001*10^16]
  R3: [1*10^17,        1.0001*10^17      ]   R4: [1*10^18, 1.0001*10^18]
  R5: [4.003*10^18,    4.0031*10^18      ]   (R5 also extends the verified range)

## Pre-registered predictions (B = 1000 unless stated)
P1. Mean rate: for each run, s' = 1/(mean survivor rank + 1) equals
    (2/ln N) / PROD_{3<=ell<=997}(1 - 1/ell), within a window-level 95% CI.
    Numerically (exact product 0.161931): 10^15: 0.3576 | 10^16: 0.3353
    | 10^17: 0.3156 | 10^18: 0.2981 | 4.003x10^18: 0.2884.
    (The B=1000 local factors measured on 2026-11-19 data were 0.9676-0.9993;
    with 10,000 patches instead of 11, lambda is predicted to lie within
    [0.995, 1.005] of 1.)
P2. Shape: survivor ranks are geometric at every scale (chi-square, window-level
    uncertainty, alpha = 0.01).
P3. Extremes (the counterexample-shaped statistic): with n = 10^6 samples and the
    run's own measured s', the global maximum survivor rank K_max satisfies
    P(K_max >= k) ~ n(1-s')^k. Numerically, median K_max ~ ln(2n)/(-ln(1-s')):
    10^15: ~33 | 10^17: ~38 | 4.003x10^18: ~43.
    REJECTION BAND: K_max above the 99.9% quantile k* = ln(1000 n)/(-ln(1-s'))
    (10^15: ~47 | 10^17: ~55 | 4.003x10^18: ~61) is a registered anomaly and will be
    investigated and reported regardless of cause. K_max below the 0.1% quantile
    likewise (too-thin tails are also a violation of H).
P4. Per-window maxima follow the same law with n = 100 (median ~ ln(200)/(-ln(1-s'))),
    and the 10,000 per-window maxima collectively fit the implied Gumbel-type
    distribution (KS on discrete-corrected quantiles, alpha = 0.01).
P5. B-scaling collapse reproduces on the new data: s'(B) proportional to the exact
    rough density for B in {100, 316, 1000, 3162, 10000}, residuals < 1%.

## Analysis
analyze_extremes.py (archived alongside this file) runs all five tests unmodified on
the five minp CSVs. Window-level block bootstrap for every uncertainty. All outcomes
reported, pass or fail.
