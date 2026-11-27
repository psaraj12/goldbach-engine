# goldbach-gaps

Work for the follow-up to the starter paper: Goldbach partitions N = p + q where
the small prime p is a *g-starter*, that is, p + g is also prime, for a fixed even gap g.

## The idea in plain words

The starter paper showed that twin starters (gap 2) and cousin starters (gap 4)
each fail on a whole third of the even numbers, for a reason that comes down to
dividing by 3. Gaps that are multiples of 6 have no such problem. Here we ask
three questions:

1. **Which gaps work?** Our conjecture: a gap g leaves only finitely many even
   numbers without a g-starter partition exactly when 6 divides g. The "only if"
   half is proved; the "if" half is open.
2. **How large must the starter be?** For each gap we measure the largest
   "smallest starter" needed, and how it grows with the size and the number of
   even numbers checked.
3. **Which gap is most efficient?** Larger gaps that contain more small primes
   (30, 210, 2310, ...) have more starters, but each one is less likely to be
   prime-paired at small sizes. The best gap depends on how far the starters reach.

## Results

See [RESULTS.md](RESULTS.md): census of 100 gaps, and two pre-registered window
experiments (Round 1: two of four criteria failed; Round 2: all five hypotheses supported).

## Contents

| Path | What it is |
|---|---|
| `engine/goldbach_v6_4_gap.cpp` | The verifier. Default build: twin/cousin starters. `-DGAP=g` (6 \| g): g-starters. Appends one summary row per finished run to `goldbach_gap<g>_minpstats.csv`. |
| `census/` | Exhaustive census of all 100 gaps g = 6, 12, ..., 600 over every even N up to 10^7 (`gap_census_1e7.csv`), and the largest minimal starter at X = 10^3 ... 10^8 (`gap_growth_1e8.csv`), with the scripts that produce and fit them. |
| `prereg/` | **Pre-registered** window experiment: frozen predictions of four candidate growth laws (`predictions.csv`), the hypotheses and pass/fail criteria (`PREREG.md`), the script that wrote them (`prereg_windows.py`), and the script that runs every window (`run_grid.sh`). |
| `runs/` | Round 1 results: one `minpstats` row per window (420 windows) and the engine log of each run. |
| `prereg_r2/` | **Round 2 pre-registration**: laws with the worst-N residue correction F(g), calibrated on Round 1 and frozen (`PREREG_R2.md`, `predictions_r2.csv`), the script (`prereg_r2.py`, which also scores the mechanism tests H4 and H5) and the runner (`run_grid_r2.sh`). |
| `runs_r2/` | Round 2 results: one `minpstats` row per window (560 windows) and the engine log of each run. |

**Not here:** the campaign source of the arXiv starter paper. It lives unchanged in
`../goldbach-starters/` and must stay byte-identical there, because the paper
states that it rebuilds the campaign binary bit for bit.

## Pre-registration

Round 1: `prereg/PREREG.md` and `prereg/predictions.csv` were committed and pushed
(commit `667f1ad`, tag `prereg-windows-v1`) before any registered window was run.
Round 2: `prereg_r2/PREREG_R2.md` and `prereg_r2/predictions_r2.csv` are committed and
pushed (tag `prereg-windows-r2`) before any Round-2 window is run. The predictions'
SHA-256 is recorded in `PREREG.md`. The model constants were calibrated on the
census only and are never refitted on window data; any later model is reported
as exploratory.

## Reproduce

From this directory:

```bash
# 1. regenerate the predictions (must reproduce the SHA-256 in PREREG.md)
python3 prereg/prereg_windows.py predict --census census/gap_growth_1e8.csv --threads 12 --out prereg

# 2. build one engine per gap (commands also listed at the top of run_grid.sh)
for g in 6 30 210 2310 30030 510510; do
  g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX2 -DGAP=$g engine/goldbach_v6_4_gap.cpp -o engine/gb_gap$g
done

# 3. run every registered window (restartable; finished windows are skipped)
ENGINE_DIR=engine THREADS=12 prereg/run_grid.sh

# 4. score the results against the frozen predictions
python3 prereg/prereg_windows.py compare --pred prereg/predictions.csv --stats-dir runs
```

Round 2, also from this directory (after its registration commit is pushed):

```bash
for g in 6 12 42 66 330 390 2310 4620; do
  g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX2 -DGAP=$g engine/goldbach_v6_4_gap.cpp -o engine/gb_gap$g
done
ENGINE_DIR=engine THREADS=12 prereg_r2/run_grid_r2.sh
python3 prereg_r2/prereg_r2.py compare --pred prereg_r2/predictions_r2.csv --stats-dir runs_r2
```

The census can be recomputed with `python3 census/gap_census.py 10000000` and
`python3 census/gap_growth.py 100000000` (about 1 minute and 20 minutes).
