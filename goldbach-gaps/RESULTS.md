# Results: g-starter Goldbach partitions and their hardest cases

## In plain words

Every even number we tested can be written as a prime plus a smaller prime p, where p
also has a prime partner a fixed distance g above it (p and p + g both prime). We call such
a p a *g-starter*.

1. **Which gaps work.** When g is a multiple of 6, only a handful of small even numbers
   (at most up to 182) lack such a partition in our census. When g is not a multiple of 6, a
   whole third of the even numbers fail, for a simple reason involving division by 3.
2. **How large the starter must be.** For each even number we take the *smallest* g-starter
   that works, and for each stretch of numbers we record the largest of these. That
   largest value grows slowly and predictably with the size of the numbers and with how
   many of them are checked.
3. **Which numbers are hardest.** The even numbers that need the largest starters are not
   random. Almost all of them sit in the residue classes, modulo small primes, that leave
   the fewest prime partners. A correction built from this, with no fitted parameters,
   explains why some gaps need larger starters than others.

Points 2 and 3 were tested in two **pre-registered** rounds: the predictions were frozen,
committed and pushed to GitHub before any of the tested windows was run. Round 1 failed
two of its four tests; its failure pointed to the residue mechanism; Round 2 tested that
mechanism on new gaps and new heights, and all five of its hypotheses were supported.

Nothing here is a proof. These are measured regularities, predicted in advance and then
checked, together with the heuristic (Hardy-Littlewood) reasoning that explains them.

---

## 1. Definitions

| Symbol | Meaning |
|---|---|
| g-starter | prime p with p + g prime (g even) |
| E(g) | even N >= 6 with no partition N = p + q, q prime, p a g-starter |
| window | W consecutive even N starting at height N |
| M | largest, over the window, of the smallest g-starter partner of each N |
| hardest N | the N in the window that attains M |
| K_g(y) | exact number of g-starters <= y |
| S(g) | Hardy-Littlewood density of g-starters relative to twin starters: product over odd primes l dividing g of (l-1)/(l-2) |
| F(g) | worst-case residue factor: product over primes l >= 5 not dividing g of (1 - 1/(l-2)^2) |
| H_g(N) | conditional factor of N: product over odd primes l <= 97 of (1 - b/l)/((1 - a/l)(1 - 1/l)), a = \|{0, -g} mod l\|, b = \|{0, -g, N} mod l\| |

"Unfavourable at l" means N is not 0 mod l, and also not -g mod l when l does not divide g.

## 2. Census of 100 gaps (exhaustive)

Every gap g = 6, 12, ..., 600, every even N from 6 to 10^7 (about 5 x 10^8 gap-N pairs),
all starters below N, no anchor limit (`census/`).

| Measure | Result |
|---|---|
| Exceptions per gap | 1 to 12; 29 gaps have N = 6 as the only exception |
| Largest exception, any gap | **182** (g = 486); none in (2 x 10^6, 10^7] for any gap |
| Starter counts vs Hardy-Littlewood | measured / predicted between 0.991 and 1.001 |
| Largest minimal starter | 1,549 (g = 30) to 4,451 (g = 222, N = 6,249,262) |

Exceptions and witnesses were re-checked by an independent script (sympy).

**Conjecture C.** For even g, E(g) is finite if and only if 6 | g. The direction
"6 does not divide g implies E(g) infinite" is proved (the blocked class is N = g mod 3);
the other direction is open.

## 3. Pre-registration record

| | Round 1 | Round 2 |
|---|---|---|
| Registration commit | `667f1ad`, 2026-11-27 23:22:29 +0530 | `0ea06a6`, 2026-11-27 23:56:08 +0530 |
| Tag | `prereg-windows-v1` | `prereg-windows-r2` |
| Frozen predictions | `prereg/predictions.csv`, SHA-256 `2c7c986d04775ac4f01fbce29fa37889d8553387f0576297ca8fef5fd3c27f5b` | `prereg_r2/predictions_r2.csv`, SHA-256 `f4df2bff5d8d2197fbdb74ba2a7f4ab2e9f849017499acddc81370d423037b3b` |
| Calibration data | census only (`census/gap_growth_1e8.csv`) | Round-1 windows only |
| Result files dated | 23:28 to 23:41 (after registration) | 23:58 to 00:18 (after registration) |
| Windows | 420 (6 gaps x 5 heights x 14) | 560 (8 gaps x 5 heights x 14) |
| Even N examined | 744,000,000 | 992,000,000 |

Every window row was checked: registered range, correct engine version, not resumed, a
partner found for all W numbers, anchor limit 10^6, hardest N inside the window.
Engine: `engine/goldbach_v6_4_gap.cpp` with `-DGAP=g`, verified against an independent
brute force at 10^12 and at 8 x 10^18.

## 4. Round 1: four candidate laws

Laws of the form size(M) = C x count, C fitted on the census:
size K = K_g(M) or S(g) Li2(M); count NW = ln N ln W or NN = (ln N)^2. Primary: K-NW.

| Law | actual/pred | mean \|log error\| | in 80% band |
|---|---|---|---|
| K-NW (primary) | 1.187 | 0.185 | 68.1% |
| L-NW | 1.268 | 0.295 | 52.9% |
| K-NN | 0.401 | 0.916 | 5.0% |
| L-NN | 0.370 | 1.014 | 4.8% |

| Hypothesis | Verdict |
|---|---|
| H1: window size matters (each NW law beats its NN partner) | **SUPPORTED** |
| H2: the K-NW drift continues with height (1.087 at 10^9, 1.237 at 4 x 10^18) | **SUPPORTED** |
| Accuracy (ratio 0.85-1.15 and 65-95% in band) | **no law adequate** |
| Gap uniformity of the best law (each gap 0.80-1.25) | **FAIL** (g = 6 at 1.394) |

The heuristic has the right shape (both ln N and ln W matter) but under-predicts the
largest starters, most of all for g = 6.

## 5. Exploratory analysis after Round 1 (not a test)

**Gap pattern.** If the hardest N tends to be unfavourable modulo each small prime, a gap
pays a factor 1 - 1/(l-2)^2 for every small prime l that does not divide it. Relative to
g = 510510, with no fitted parameter:

| g | observed relative under-prediction | worst-case penalty prediction |
|---|---|---|
| 6 | 1.245 | 1.202 |
| 30 | 1.088 | 1.068 |
| 210 | 1.034 | 1.025 |
| 2310 | 1.007 | 1.013 |
| 30030 | 1.010 | 1.004 |

**Where the hardest N sit** (420 Round-1 windows): unfavourable at 3, 5, 7, 11 and 13
simultaneously in 88.8% of windows, against 31.6% for random even N. These observations
motivated Round 2, which tests them on new data.

## 6. Round 2: the residue mechanism, tested

New gaps (12, 42, 66, 330, 390, 4620, plus 6 and 2310), new heights
(3 x 10^9 to 8 x 10^18). Laws fitted once on Round 1:

| Law | form | actual/pred | mean \|log error\| | in 80% band | adequate |
|---|---|---|---|---|---|
| **R2-F (primary)** | K_g(M) = 0.6459 ln N ln W / F(g) | **0.983** | **0.103** | **80.2%** | **yes** |
| R2-FH | K_g(M) = 0.4320 (ln N)^1.1157 ln W / F(g) | 0.979 | 0.100 | 79.6% | yes |
| K-NW (baseline) | K_g(M) = 0.6883 ln N ln W | 1.059 | 0.129 | 74.1% | yes |

| Hypothesis | Result | Verdict |
|---|---|---|
| H1: F(g) improves prediction | mean \|log error\| 0.103 vs 0.129 | **SUPPORTED** |
| H2: height law holds (R2-FH, top / bottom height within 0.90-1.10) | 0.969 | **SUPPORTED** |
| H3: same prime divisors behave alike (within 0.90-1.10) | 12 vs 6: 1.067; 4620 vs 2310: 0.995 | **SUPPORTED** |
| H4: hardest N unfavourable at l = 3, 5, 7, 11, 13 (p < 0.001 each) | 560, 560, 551, 542, 529 of 560 windows (z = +16.7 to +6.0) | **SUPPORTED** |
| H5: hardest N have small H_g(N) (mean rank < 0.40, sign test p < 0.001) | mean rank 0.129; 559 of 560 windows below the median | **SUPPORTED** |

R2-F by gap: 6: 0.934, 12: 0.996, 42: 0.990, 66: 0.984, 330: 0.985, 390: 0.976,
2310: 1.001, 4620: 0.996. By height: 0.934 (3 x 10^9), 0.970, 1.004, 1.000, 1.006 (8 x 10^18).

## 7. Limitations and open points

1. **g = 6 remains the outlier** (0.934 in Round 2; it was also the outlier in Round 1).
   Gaps 6 and 12 share their prime divisors but not their smallest starters (5, 7, 11, 13
   against 5, 7, 11, 17, 19), so finite-range starter structure, possibly correlations
   between starters, is not captured.
2. **Mild drift at the lowest heights.** A power law in ln N and an affine law
   (alpha ln N - beta, the form Oliveira e Silva et al. found for classical Goldbach) fitted
   to the same data differ by less than 0.1% below 2^64, so their long-range behaviour cannot
   be separated by 64-bit computation.
3. **Extreme-value scatter.** Replicate windows of equal gap, height and size scatter by
   about 0.14 in ln M; individual windows are predicted to about 10%, the overall level to about 2%.
4. **Heuristic, not proof.** The laws rest on the Hardy-Littlewood framework; Conjecture C
   (sufficiency) and the growth law remain conjectural.

## 8. Related work

- Z. Juhasz and M. Bartalos, arXiv:2510.21870 (2025): ranks minimal primes in
  partitions n = m1 p + m2 q by a density factor R(m1, m2); the same factor ranks averages
  and maxima; absolute sizes left open. Our S(g) plays the role of their R; the worst-case
  factor F(g) and the location of the hardest N are not part of their framework.
- A. Granville, J. van de Lune and H. J. J. te Riele (1989): p(n) << log^2 n log log n for
  classical Goldbach; the distribution of minimal primes governed by pair factors D(q - r).
  Our heuristic with all primes as anchors reproduces this growth form; for g-starters it
  predicts one extra log log factor.
- T. Oliveira e Silva, S. Herzog and S. Pardi, Math. Comp. 83 (2014): records of p(n) to
  4 x 10^18 (maximum 9,781); affine decay law 0.755 log x - 4.19; note that extreme values
  may behave differently from typical ones.

## 9. Reproduce

From `goldbach-gaps/`:

```bash
# census
python3 census/gap_census.py 10000000
python3 census/gap_growth.py 100000000

# engines (one per gap)
for g in 6 12 30 42 66 210 330 390 2310 4620 30030 510510; do
  g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX2 -DGAP=$g engine/goldbach_v6_4_gap.cpp -o engine/gb_gap$g
done

# Round 1
ENGINE_DIR=engine THREADS=12 prereg/run_grid.sh
python3 prereg/prereg_windows.py compare --pred prereg/predictions.csv --stats-dir runs

# Round 2
ENGINE_DIR=engine THREADS=12 prereg_r2/run_grid_r2.sh
python3 prereg_r2/prereg_r2.py compare --pred prereg_r2/predictions_r2.csv --stats-dir runs_r2
```
