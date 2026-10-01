# goldbach-engine

Verification of the even Goldbach conjecture with **cached partitions** (QHot):
a sieve-free verifier that keeps a small ring of recently successful large primes
*q* and certifies up to 64 consecutive even numbers per 64-bit operation.

**Result.** Every even integer in [4·10^18, 4.001·10^18] is a sum of two primes
(5·10^14 numbers, zero exceptions). Together with Oliveira e Silva, Herzog and
Pardi (Math. Comp. 83, 2014), the even Goldbach conjecture holds up to 4.001·10^18.
The whole range was recomputed independently with a second build and different
parameters; all minimal-partition witnesses agree.

- Paper: *to be added*
- Data archive (checkpoints, logs, provenance): Zenodo, doi:10.5281/zenodo.*to be assigned*

## Contents

| Path | What it is |
| --- | --- |
| `src/goldbach_v6_2.cpp` | **v6.2, the current verifier**: v6.1 plus words carried between batches in the prefix (AVX-512, AVX2 and NEON; halves the prefix's loads); K=56 with AVX-512 |
| `src/goldbach_v6_1.cpp` | v6.1: batched QHot, sorted ring, fixed prefix (K=40 with AVX-512, else 32), optional AVX-512 / AVX2 / NEON prefix, fully-covered-batch shortcut |
| `src/goldbach_v6_1_portable.cpp` | v6.1 with the vector prefix chosen at run time (portable binaries; carried words need per-thread state and are not in the portable build) |
| `src/goldbach_v6_simd.cpp` | v6 (K=32 default, no shortcut): the build used for the independent recomputation and the 4.001→4.002·10^18 run |
| `src/goldbach_v5_batch.cpp` | v5, the build that produced the 4·10^18 → 4.001·10^18 campaign |
| `src/goldbach_v4.cpp` | v4, per-number QHot (reference implementation for equivalence tests) |
| `baselines/goldbach_v4_window.cpp` | double-sieve-style baseline (primesieve window + shift-OR), untuned |
| `experiments/` | measured variants kept for the record (see `experiments/README.md`) |
| `scripts/` | build, certificate and witness verification, run comparison |
| `docs/RESULTS.md` | campaign and performance tables |
| `docs/PROVENANCE.md` | checksums tying results to exact sources and binaries |

## Build (Linux, GCC ≥ 11)

Recommended: build on the machine that will run, with profile-guided optimization
and all checks, in one step:

```bash
sudo apt install -y g++ python3-sympy
bash scripts/build_v6.sh          # builds src/goldbach_v6_2.cpp -> build/goldbach_v6 and build/provenance_v6.txt
```

Manual build (no PGO), choosing the vector prefix explicitly:

```bash
g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX512 src/goldbach_v6_2.cpp -o goldbach_v6_2
#   -DSIMD_AVX2 (x86 without AVX-512), -DSIMD_NEON (ARM64, use -mcpu=native), or none (scalar)
```

Portable binary (one file for any x86-64 Linux machine since about 2013; picks
AVX-512, AVX2 or scalar at start-up, a few percent slower than a native build):

```bash
g++ -O3 -march=x86-64-v3 -fopenmp -std=c++17 -static src/goldbach_v6_1_portable.cpp -o goldbach_v6_1_linux_x86-64-v3
```

macOS on Apple Silicon (Homebrew GCC; `-mcpu=apple-m1` runs on every M-series Mac):

```bash
g++-16 -O3 -mcpu=apple-m1 -fopenmp -std=c++17 src/goldbach_v6_1_portable.cpp -o goldbach_v6_1_macos_arm64
```

## Run

```
goldbach_v6 START END BLOCK_BITS THREADS ANCHOR_LIMIT SAMPLE_LIMIT RING
goldbach_v6 --resume 0 0 BLOCK_BITS THREADS ANCHOR_LIMIT SAMPLE_LIMIT RING
```

- `START`, `END`: even integers, inclusive range.
- Recommended settings: `BLOCK_BITS=24`, `ANCHOR_LIMIT=40000000`, `RING=512` (power of two),
  `SAMPLE_LIMIT=250000` (total witness rows kept for the whole run).
- `THREADS`: use all hardware threads (`$(nproc)`).

Example: re-verify 10^11 even numbers just above 4·10^18

```bash
./goldbach_v6 4000000000000000000 4000000199999999998 24 $(nproc) 40000000 1000 512
```

The program writes an atomic checkpoint after every 10^11 even numbers, with
counters, a configuration record and witness samples (`N,p,q,source` rows; one
QHOT and one COLD witness per 10^10 integers). COLD witnesses are minimal
partitions. Any number without a partition from the anchor primes is escalated
and recorded in the miss file; the exit status is non-zero if misses occur.

## Verify

```bash
# every witness in a checkpoint (SymPy, deterministic below 2^64); COLD rows also checked for minimality
python3 scripts/verify_witnesses.py goldbach_v62_checkpoint.csv --minimal

# COLD witnesses of two runs over the same range must be identical
python3 scripts/compare_witnesses.py reference_checkpoint.csv rerun_checkpoint.csv

# full certificate dump for a small range (debug build), checked row by row
g++ -O3 -march=native -fopenmp -std=c++17 -DDUMP_ALL src/goldbach_v6_2.cpp -o goldbach_v6_dump
S=4000000000000000000; E=$((S + 2*1048576 - 2))
./goldbach_v6_dump $S $E 24 $(nproc) 40000000 1000 512 > /dev/null
python3 scripts/verify_dump_parallel.py goldbach_v62_dump.csv $S $E 40000000

# bit-exact check of a vector prefix against the scalar computation
g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX512 -DSIMD_CHECK src/goldbach_v6_2.cpp -o check
./check $S $((S + 2*33554432 - 2)) 24 $(nproc) 40000000 1000 512 | grep "SIMD checks"
```

## Method in one paragraph

For 64 consecutive even numbers N0, N0+2, …, N0+126 and a cached prime q, the
partners N − q are 64 consecutive odd numbers, so one 64-bit slice of a bitset of
small primes (the anchors, up to 4·10^7) tells which of the 64 numbers q certifies.
Slices from the ring of cached q are OR-ed until the batch is covered; numbers left
over are handled by an ascending scan over the anchors with deterministic
Miller–Rabin, and any new q found there joins the ring. Fewer than 1 in 100,000
numbers reach that path. The first 56 probes (with AVX-512 or NEON; 32 otherwise)
are a fixed, unrolled, optionally vectorized prefix, following the per-processor
tuned inner loop of Oliveira e Silva et al. (Algorithm 1.4), and batches fully
covered by the ring skip the per-number loop entirely. With AVX-512 or AVX2, each
entry's second bitset word is kept for the next batch, where it is the first word
needed, halving the prefix's loads.

## Acknowledgments

The fully-covered-batch shortcut (v6.1) and the reuse of prefix addressing between
<<<<<<< HEAD
batches, developed into carried words (v6.2), were suggested in code reviews by ChatGPT. The fixed-length prefix follows the tuned inner loop of Oliveira e Silva,
Herzog and Pardi (2014).
=======
batches, developed into carried words (v6.2), were suggested in code reviews by
ChatGPT. The code was developed with Claude. The fixed-length prefix follows the
tuned inner loop of Oliveira e Silva, Herzog and Pardi (2014).
>>>>>>> 2b9a58f (README: update method paragraph and acknowledgments)

## Citation

See `CITATION.cff`. Please cite the paper and the Zenodo archive.

## License

Code: MIT (see `LICENSE`). Data in the Zenodo archive: CC-BY 4.0.
