# goldbach-engine

A sieve-free verifier for the even Goldbach conjecture. It checks about
**4·10^10 even numbers per second** on a single 24-OCPU cloud VM, which is about
1.5 CPU cycles per even number per core.

**Result.** Every even integer in [4·10^18, 4.003·10^18] is a sum of two primes:
1.5·10^15 numbers, with zero exceptions. Combined with Oliveira e Silva, Herzog and
Pardi (Math. Comp. 83, 2014), this means the even Goldbach conjecture holds up to
4.003·10^18. The latest 5·10^14 numbers (4.002 → 4.003·10^18) took 3 h 21 min on
one AMD EPYC 9J14 VM.

**What is different.** Earlier large verifications find the large prime *q* in
N = p + q by sieving the numbers just below N. goldbach-engine never sieves large
numbers. It keeps a small ring of recently successful primes *q* (the "Q-Hot cache").
For 64 consecutive even numbers, one 64-bit word of a small-prime bitset shows which
of them a cached *q* certifies. More than 99.999% of numbers are certified this way.
The rest fall back to an ascending scan over small primes *p*, with deterministic
Miller–Rabin used to prove N − p prime.

- Method paper: S. A. R. Parthibanathan, *Q-Hot Cache: A High-Throughput Method for
  Empirical Verification of the Even Goldbach Conjecture*, Zenodo preprint, April 2026,
  [doi:10.5281/zenodo.19884541](https://doi.org/10.5281/zenodo.19884541).
- Data (checkpoints, logs and provenance for every run):
  [doi:10.5281/zenodo.23082138](https://doi.org/10.5281/zenodo.23082138).

## Verification status

| Range | Even integers | Primary run | Independent rerun | Exceptions |
| --- | --- | --- | --- | --- |
| [4.000, 4.001]·10^18 | 5.0·10^14 | v5 | v6 and v6.2 (and earlier, campaign 1) | 0 |
| [4.001, 4.002]·10^18 | 5.0·10^14 | v6 (AVX-512, K=40) | v5 | 0 |
| [4.002, 4.003]·10^18 | 5.0·10^14 | v6.2 (AVX-512, K=56) | v5 | 0 |

Every part of the range has been verified by two different programs, and each run
has its source, binary, log and checkpoint on record. The minimal-partition (COLD)
witnesses of each pair of overlapping runs were compared: all 300,000 are
identical. All 1,200,000 recorded witnesses from these runs pass an independent
SymPy check. For per-run hardware, timing and hit rates, see
[`docs/RESULTS.md`](docs/RESULTS.md). For checksums that tie each result to its
exact source, binary and PGO profile, see [`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## Why the result is a proof for the range

A number N is counted as verified only when the program has exhibited N = p + q with
both p and q proven prime:

- **p** comes from a table of primes up to 4·10^7, built with an ordinary sieve of
  Eratosthenes.
- **q** is proven prime by Miller–Rabin with the 7 bases {2, 325, 9375, 28178,
  450775, 9780504, 1795265022} (Sinclair). This test is deterministic for all
  n < 2^64, and it has been checked against Feitsma's database of base-2
  pseudoprimes. Bases a ≥ n are skipped rather than reduced mod n. Every n at or
  below the largest base is either removed by trial division or decided by the
  smaller bases. A *q* enters the Q-Hot ring only after passing this test, so a
  ring hit is a complete certificate.
- Any N that the anchor primes do not cover is escalated, recorded in a miss file,
  and makes the program exit with a non-zero status. All campaigns ended with zero
  misses.

Correctness therefore rests on the program and its execution, as it does in every
large computational verification. The evidence for the computation is:

- independent reruns with different programs, whose COLD witnesses must agree exactly
- bit-exact checks of every vector (SIMD) path against the scalar computation
- a full certificate dump for small ranges, checked row by row
- one QHOT and one COLD witness sampled per 10^10 integers in every campaign,
  re-checked with SymPy

The sampled witnesses are a spot check of the run. They are not a certificate for
every number in it.

## Quick start

You need Linux (or WSL on Windows) with GCC 11 or newer, and Python 3 with SymPy.
For macOS, see [Build](#build).

**1. Get the code and build it.** This takes about 5 minutes and runs all
correctness checks.

```bash
sudo apt install -y git g++ python3-sympy
git clone https://github.com/psaraj12/goldbach-engine.git
cd goldbach-engine
bash scripts/build_v6.sh
```

The script picks the best vector instructions your CPU supports (AVX-512, AVX2 or
scalar). It builds a profile-guided binary at `build/goldbach_v6`, and stops with a
message if any check fails.

**2. Verify a range.** This command checks 10^11 even numbers just above 4·10^18.
That takes a few seconds on a large server and under a minute on a recent laptop.

```bash
mkdir -p run && cd run
../build/goldbach_v6 4000000000000000000 4000000199999999998 24 $(getconf _NPROCESSORS_ONLN) 20000000 1000 512
```

At the end, you should see `Coverage: N / N` and `Misses: 0`.

**3. Check the witnesses independently.**

```bash
python3 ../scripts/verify_witnesses.py goldbach_v63_checkpoint.csv --minimal
```

This re-checks every recorded witness N = p + q with SymPy. It also confirms that
each COLD witness uses the smallest possible p.

## Run

```
goldbach_v6 START END BLOCK_BITS THREADS ANCHOR_LIMIT SAMPLE_LIMIT RING
goldbach_v6 --resume 0 0 BLOCK_BITS THREADS ANCHOR_LIMIT SAMPLE_LIMIT RING
```

| Argument | Meaning | Use |
| --- | --- | --- |
| `START`, `END` | even integers, inclusive | |
| `BLOCK_BITS` | work-block size (log2) | `24` |
| `THREADS` | worker threads | `$(getconf _NPROCESSORS_ONLN)` |
| `ANCHOR_LIMIT` | largest small prime p | `20000000` (v6.3; the campaigns used `40000000`) |
| `SAMPLE_LIMIT` | cap on witness rows kept for the whole run | `1000` for tests, `250000` for campaigns |
| `RING` | Q-Hot ring size (power of two) | `512` |

The binary's PGO profile is trained on these settings, so only `START`, `END`,
`THREADS` and `SAMPLE_LIMIT` normally change. The program writes an atomic
checkpoint after every 10^11 even numbers. Each checkpoint contains counters, a
configuration record and witness rows (`N,p,q,source`). To resume an interrupted
run, use `--resume` from the same folder. Run only one instance per folder.

## Build

**Linux, recommended.** Build on the machine that will do the run. This uses
profile-guided optimization and runs all checks:

```bash
bash scripts/build_v6.sh      # src/goldbach_v6_3.cpp -> build/goldbach_v6, build/provenance_v6.txt
```

**Manual build** (no PGO), with the vector prefix chosen explicitly:

```bash
g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX512 src/goldbach_v6_3.cpp -o goldbach_v6_3
#   -DSIMD_AVX2 (x86 without AVX-512), -DSIMD_NEON (ARM64, use -mcpu=native), or none (scalar)
```

Leave out `-DCARRYIDX` in a manual build. The carried index pays off only together with
PGO: with PGO it adds about 5% on AVX2 and 3–13% on AVX-512, but without PGO it cost
about 4% on AVX2. `scripts/build_v6.sh` enables it, because that script always uses PGO.

**Portable binary.** This builds one file that runs on any x86-64 Linux machine
from about 2013 on. It chooses AVX-512, AVX2 or scalar at start-up and is a few
percent slower than a native build.

```bash
g++ -O3 -march=x86-64-v3 -fopenmp -std=c++17 -static src/goldbach_v6_1_portable.cpp -o goldbach_v6_1_linux_x86-64-v3
```

**macOS (Apple Silicon).** Use the macOS build script, which does PGO and runs all
checks:

```bash
brew install gcc
python3 -m venv .venv && source .venv/bin/activate
pip install sympy
bash scripts/build_v6_mac.sh
```

## Checking a build or a run

```bash
# every witness in a checkpoint; COLD rows are also checked for minimality
python3 scripts/verify_witnesses.py goldbach_v63_checkpoint.csv --minimal

# COLD witnesses of two runs over the same range must be identical
python3 scripts/compare_witnesses.py reference_checkpoint.csv rerun_checkpoint.csv

# full certificate dump for a small range (debug build), checked row by row
g++ -O3 -march=native -fopenmp -std=c++17 -DDUMP_ALL src/goldbach_v6_3.cpp -o goldbach_v6_dump
S=4000000000000000000; E=$((S + 2*1048576 - 2))
./goldbach_v6_dump $S $E 24 $(getconf _NPROCESSORS_ONLN) 20000000 1000 512 > /dev/null
python3 scripts/verify_dump_parallel.py goldbach_v63_dump.csv $S $E 20000000

# bit-exact check of the vector prefix against the scalar computation
g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX512 -DCARRYIDX -DSIMD_CHECK src/goldbach_v6_3.cpp -o check
./check $S $((S + 2*33554432 - 2)) 24 $(getconf _NPROCESSORS_ONLN) 20000000 1000 512 | grep "SIMD checks"
```

## Method in more detail

Take 64 consecutive even numbers N0, N0+2, …, N0+126 and a cached prime q. The
partners N − q are then 64 consecutive odd numbers. So a single 64-bit slice of the
small-prime bitset shows which of the 64 numbers q certifies. Slices from the
cached q are OR-ed together until every number in the batch is covered. If a whole
batch is covered by the ring, the per-number loop is skipped. Any numbers left over
go to the ascending anchor scan with Miller–Rabin, and each new q found there joins
the ring.

The first 56 probes (AVX-512 or NEON; 32 otherwise) form a fixed, unrolled,
vectorized prefix. This follows the tuned inner loop of Oliveira e Silva et al.
(Algorithm 1.4). With AVX-512 or AVX2, the second bitset word of each ring entry is
carried over to the next batch, where it becomes the first word needed. This halves
the prefix's memory loads. v6.3 also carries the gather index and shifts (with PGO),
and drops work from the cold path that the ring walk makes redundant: a number that
reaches the cold path cannot have its q in the ring, and the trial division by 2–17
already done is not repeated before Miller–Rabin.

Throughput on one 12-OCPU EPYC 9J14 VM, for the same work:

| Version | Change | M evens/s |
| --- | --- | --- |
| v4 | Q-Hot per number | 915 |
| v5 | batched Q-Hot | 5,494 |
| v6 | sorted ring, PGO, AVX-512 prefix | 9,618 |
| v6.1 | fully-covered-batch shortcut | 14,825 |
| v6.2 | carried words | 18,641 |

On 24 OCPUs, v6.2 with K=56 reaches 41,526 M/s. v6.3 gives identical output and, built
with PGO, is faster than v6.2 on x86 and level with it on ARM: one thread, +6.9% on an
Alder Lake laptop (AVX2) and +0.3% on an Apple M4 (NEON). Results for more machines are
in [`docs/RESULTS.md`](docs/RESULTS.md).

## Repository layout

| Path | What it is |
| --- | --- |
| `src/goldbach_v6_3.cpp` | **current verifier (v6.3)**: v6.2 plus a lean cold path and, with PGO, the carried index; identical output |
| `src/goldbach_v6_2.cpp` | v6.2; produced the 4.002→4.003·10^18 run and the recomputation of segment 2 |
| `src/archive/goldbach_v6_2_11525e8d.cpp` | the exact revision used for those runs (rebuilds the campaign binary bit-identically) |
| `scripts/` | build scripts, witness and certificate verification, run comparison |
| `docs/RESULTS.md` | campaign and performance tables |
| `docs/PROVENANCE.md` | SHA-256 checksums tying each result to its source, binary and profile |
| `docs/provenance/` | the original provenance record of every run, and checksums of all archived data |
| `data/README.md` | layout of the Zenodo data record and how to check it |
| `baselines/` | a double-sieve-style baseline (primesieve window + shift-OR), untuned |
| `experiments/` | measured variants kept for the record (see `experiments/README.md`) |

<details>
<summary>Earlier versions (used for the campaigns and reruns)</summary>

| Path | What it is |
| --- | --- |
| `src/goldbach_v6_1.cpp` | v6.1: batched Q-Hot, sorted ring, fixed prefix (K=40 with AVX-512, else 32), fully-covered-batch shortcut |
| `src/goldbach_v6_1_portable.cpp` | v6.1 with the vector prefix chosen at run time (no carried words) |
| `src/goldbach_v6_simd.cpp` | v6: the independent recomputation and, with K=40, the 4.001→4.002·10^18 run |
| `src/goldbach_v5_batch.cpp` | v5: the 4·10^18 → 4.001·10^18 campaign |
| `src/goldbach_v4.cpp` | v4, per-number Q-Hot: the program generation of the first verification of 4·10^18 → 4.001·10^18, and the reference implementation for equivalence tests |

</details>

## Author

Santhiagu Antony Raj Parthibanathan ([ORCID 0009-0002-0355-4191](https://orcid.org/0009-0002-0355-4191)).
Idea, design decisions, computations and verification.

## Acknowledgments

Implementation, benchmarking and checking were carried out with the assistance of
Claude (Anthropic). Code reviews by ChatGPT (OpenAI) suggested the
fully-covered-batch shortcut (v6.1) and the reuse of prefix addressing between
batches, which was developed into carried words (v6.2). The fixed-length prefix
follows the tuned inner loop of Oliveira e Silva, Herzog and Pardi (2014).

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite both the paper
(doi:10.5281/zenodo.19884541) and the data record (doi:10.5281/zenodo.23082138).

## License

Code: MIT (see `LICENSE`). Data in the Zenodo archive: CC-BY 4.0.
