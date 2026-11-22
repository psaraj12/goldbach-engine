# goldbach-starters

Code, data and provenance for the preprint

> **Goldbach Partitions with a Twin, Cousin or Sexy Starter: a Mod-3 Obstruction,
> a Density Law, and Verification to 10¹³**
> Santhiagu Antony Raj Parthibanathan (ORCID [0009-0002-0355-4191](https://orcid.org/0009-0002-0355-4191))

The paper is in [`paper/`](paper/) (LaTeX source and PDF).

## The idea in plain words

Goldbach's conjecture says every even number from 4 upward is a sum of two primes.
In practice, one of the two primes can almost always be chosen very small, and that small
prime can almost always be chosen next to another prime, such as a p with p + 2 also prime.
This repository verifies two precise versions of that observation and contains everything
needed to check the verification.

* A **twin starter** is a prime p with p + 2 prime, a **cousin starter** has p + 4 prime,
  and a **sexy starter** has p + 6 prime.
* **Conjecture A:** every even N ≥ 6 is p + q with q prime and p a twin or cousin starter.
* **Conjecture B:** every even N ≥ 8 is p + q with q prime and p a sexy starter.

## Results

| Conjecture | Range (all even N) | Even N checked | Misses | Rate |
|---|---|---|---|---|
| A | [6, 10¹²] | 499,999,999,998 | 0 | 502.7 M/s |
| A | [10¹² + 2, 10¹³] | 4,500,000,000,000 | 0 | 519.5 M/s |
| B | [8, 10¹²] | 499,999,999,997 | 0 | 622.5 M/s |

Anchor limits were P = 2·10⁷ for the first campaign and P = 10⁷ for the others (every
witness has p ≤ P). All three ran on a laptop with an
Intel Core i5-1235U (12 threads, WSL2), with GCC 13.3, AVX2 and profile-guided builds.

## Repository layout

```
src/
  goldbach_v6_4_tc.cpp     the verifier used for this paper (starter anchors);
                           -DSEXY_ANCHORS selects Conjecture B, default is Conjecture A
  goldbach_v6_4.cpp        the unrestricted Q-Hot verifier it was derived from (all odd primes)
scripts/
  build_v6.sh              PGO build + bit-exact SIMD check + forced-miss check + provenance
  tc_analyze.py            analysis of a minimal-witness log (gearbox table, density law, fuel gauge)
  audit_tc.py              independent sympy audit of a Conjecture-A campaign
  audit_sx.py              independent sympy audit of a Conjecture-B campaign
  exact_needs4_bands_1e8.py  exact counts for Table 1 up to 1e8 (numpy)
  single_gap_tests_1e7.py  single-gap versions (g = 2, 4, 6, 8, 12) to 1e7
  lowside_check.py         low-side form of Conjecture A (p <= N/2) to 1e7
  propagation_steps.py     the +2 / +4 propagation steps to 1e7
  deep_probe_1e18.py       random spot probes at 1e9 .. 1e18
  crosscheck_pandey_records.py  reproduces the 51 least-witness records of Pandey
                           (arXiv:2608.02381) and checks that for N = 1 (mod 3) they equal the
                           least twin/cousin starter (Proposition 9 of the paper)
paper/
  goldbach_starters.tex / .pdf
provenance/                compiler, flags, CPU and SHA-256 of source, binary and PGO profile
results/
  conjecture_A/campaign_6_1e12/       checkpoint, misses, minimal-witness log (500,000 rows), analysis,
                                      crawl curve, audit output, transcribed console summary
  conjecture_A/campaign_1e12_1e13/    run.log, checkpoint, misses, minimal-witness log (900,000 rows),
                                      crawl curve, analysis, audit output
  conjecture_B/campaign_8_1e12/       run.log, checkpoint, misses, minimal-witness log (500,000 rows),
                                      audit output, summary of the plain-build rerun
  exploratory/             outputs of the exploratory scripts
SHA256SUMS                 checksums of every file in this directory tree
```

## Provenance chain

| File | Source SHA-256 | Binary SHA-256 | Produced |
|---|---|---|---|
| `provenance_tc_campaign_2026-11-21.txt` | `37da1af2…` | `19ad4663…` | the binary for both Conjecture A campaigns |
| `provenance_tc_rebuild_2026-11-22.txt` | `841c219e…` | `19ad4663…` | a rebuild from the final source in `src/` |
| `provenance_sx_campaign_2026-11-22.txt` | `841c219e…` | `6d572f85…` | the binary for the Conjecture B campaign |

Rebuilding Conjecture A from the published source (`841c219e…`) gives **exactly the
campaign binary** (`19ad4663…`). All changes between the two source versions are comments
or code under `#ifdef SEXY_ANCHORS`. PGO profile files (`gs.gcda`) differ between builds
because GCC stamps them per compilation; the profiled binaries are identical.

## Reproducing

```bash
# Conjecture A (twin/cousin starters)
GB_BUILD_DIR=~/gb_build_tc bash scripts/build_v6.sh src/goldbach_v6_4_tc.cpp
~/gb_build_tc/goldbach_v6 6 1000000000000 24 $(nproc) 20000000 1000 512 --minp-per-window=5000
~/gb_build_tc/goldbach_v6 1000000000002 10000000000000 24 $(nproc) 10000000 1000 512 --minp-per-window=1000

# Conjecture B (sexy starters)
GB_SEXY=1 GB_BUILD_DIR=~/gb_build_sx bash scripts/build_v6.sh src/goldbach_v6_4_tc.cpp
~/gb_build_sx/goldbach_v6 8 1000000000000 20 $(nproc) 10000000 1000 512 --minp-per-window=5000

# analysis and audit
python3 scripts/tc_analyze.py goldbach_tc_minp.csv
python3 scripts/audit_tc.py results/conjecture_A/campaign_1e12_1e13
python3 scripts/audit_sx.py results/conjecture_B/campaign_8_1e12
```

`build_v6.sh` expects `verify_dump.py` (from the Q-Hot verifier in
[psaraj12/goldbach-engine](https://github.com/psaraj12/goldbach-engine)) next to it in
`scripts/` for its certificate-dump stage. Copy or symlink it there before building. On WSL, build and
run on the native Linux filesystem (`~/…`), not under `/mnt/*`.

## Notes on the data

* Each `audit_output.txt` is produced by `scripts/audit_tc.py` or `scripts/audit_sx.py`.
  The audits recheck every sampled witness and every logged minimal witness with sympy,
  spot-check minimality, and (for Conjecture A) confirm Proposition 3 row by row: for every
  logged N = 2 (mod 3), a twin-starter partner exists exactly when N - 3 is prime.
* **Conjecture B:** the run of record is the PGO run. A plain-build rerun reported identical
  work counters, and its minimal-witness data rows are byte-identical to those included here.
* **Conjecture A, [6, 10^12]:** no `run.log` was kept for this run. Its console summary was
  transcribed; the settings are taken from the checkpoint header.
* Campaign settings (from the checkpoint headers):

  | Campaign | Block bits | P | `--minp-per-window` |
  |---|---|---|---|
  | A, [6, 10^12] | 24 | 2·10^7 | 5000 |
  | A, [10^12 + 2, 10^13] | 24 | 10^7 | 1000 |
  | B, [8, 10^12] | 20 | 10^7 | 5000 |

  All used 12 threads, sample limit 1000 and a ring of 512.

## Related work by the author

The verifier is the Q-Hot method, with its anchor table restricted to starters:

* S. A. R. Parthibanathan, *Q-Hot Cache: a high-throughput method for empirical
  verification of the even Goldbach conjecture*, Zenodo preprint, April 2026,
  doi:[10.5281/zenodo.19884541](https://doi.org/10.5281/zenodo.19884541).
* Verification of the even Goldbach conjecture to 4.003·10¹⁸ (submitted). Campaign data
  doi:[10.5281/zenodo.23082138](https://doi.org/10.5281/zenodo.23082138), production audit
  doi:[10.5281/zenodo.23278480](https://doi.org/10.5281/zenodo.23278480).

## License

Code (`src/`, `scripts/`): MIT, see `LICENSE`. Data, provenance and paper: CC BY 4.0.

## How to cite

Please cite the preprint (arXiv identifier to be added) and this repository.
