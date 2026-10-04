# Provenance

Each result is tied by SHA-256 to the exact source file, binary and (for PGO builds)
profile that produced it. The original provenance records of every run are kept
verbatim in [`provenance/`](provenance/). [`provenance/DATA_SHA256SUMS.txt`](provenance/DATA_SHA256SUMS.txt)
lists the SHA-256 of every file in the Zenodo data record
([doi:10.5281/zenodo.23082138](https://doi.org/10.5281/zenodo.23082138)):
checkpoints, miss files, logs and provenance records.

Runs 2–7 were made on Oracle Cloud VM.Standard.E5.Flex instances (AMD EPYC 9J14). Campaign 1 (run 1) has checkpoints only; see the notes.

## Sources

| Repository file | SHA-256 | Used for |
| --- | --- | --- |
| `src/goldbach_v5_batch.cpp` | `b23c765a1ac9232d7253b005eedc018b001d1a2e6f77d383236c5ca6bcc8fb0a` | v5 campaign (both segments); v5 recomputation of 4.001→4.003·10^18 |
| `src/goldbach_v6_simd.cpp` | `6a7ace5fd6a37ba820a4b814067dc4f3582076f5d5f8b41d1ac7e46f5906f0af` | v6 recomputation of segment 1; v6 extension 4.001→4.002·10^18 |
| `src/archive/goldbach_v6_2_11525e8d.cpp` | `11525e8d03ba285d9541516399dc0ddd3732b8977b5087b4d7a4ca26c06d10d2` | v6.2 extension 4.002→4.003·10^18; v6.2 recomputation of segment 2 (same binary) |
| `src/archive/goldbach_verifier_ckpt_be0158aa.cpp` | `be0158aa0b196a2737c2e9071e1c9a29fbd1a5045d19fc830ecda41c49a8c713` | campaign 1, invocations 1 and 2 and the run behind checkpoints 3–12 (see the note below) |

Each repository file above is byte-identical to the copy archived with its run.

## Runs

| # | Run | Range (even N, inclusive) | Build | Binary SHA-256 | PGO profile SHA-256 | VM |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | campaign 1 | 4·10^18 → 4.001·10^18 (both ends included) | not recorded | not recorded | n/a | not recorded |
| 2 | v5 campaign, segment 1 | 4·10^18 → 4,000,318,999,999,999,998 | GCC 13.3.0, `-O3 -march=native -fopenmp`, no PGO | `8187c1c7d0486da02e6e48adc398cb69c9bcb4065cd8a68c4cbb38b5c1a05097` | n/a | 12 OCPUs, 24 threads |
| 2 | v5 campaign, segment 2 | 4,000,319,000,000,000,000 → 4.001·10^18 − 2 | same binary as segment 1 | `8187c1c7…5097` | n/a | 12 OCPUs, 24 threads |
| 3 | v6 recomputation of segment 1 | 4·10^18 → 4,000,318,999,999,999,998 | `-DSIMD_AVX512`, K=32, PGO | `7d8f4ec3700effb6f6e8a2606028c85c82336077aee325566ba33c7d7c481efe` | `cc3890f2dfaf08b5d11fa7a985b6e10869e357f49844625d6c1ac5e543c592d7` | 24 OCPUs, 48 threads |
| 4 | v6 extension | 4.001·10^18 → 4.002·10^18 − 2 | GCC 15.2.0, `-DSIMD_AVX512 -DKFIX=40`, PGO | `9d3cc2808307ea263b6756e92bb5b599843be1d0f7de431cc4c856dc0bd3301c` | `4f7de34b64fa5f86abd0875316d6260b60cfe520642e69f492ae02bca4a0f687` | 24 OCPUs, 48 threads |
| 5 | v6.2 extension | 4.002·10^18 → 4.003·10^18 − 2 | GCC 15.2.0, `-DSIMD_AVX512 -DKFIX=56`, PGO | `bfea4378691c270f0d12e1a4d65b1f1d059c74ce370ee0e7da723bc2b017a370` | `cd2d1cd06a0c69d4d8e239637889e3e5fa2caf7520f83f8d9fec9af194bfe683` | 24 OCPUs, 48 threads |
| 6 | v5 recomputation of the extensions | 4.001·10^18 → 4.003·10^18 − 2 | GCC 13.4.0, `-O3 -march=native -fopenmp`, no PGO | `216f17eb9ef8cf5f21bf31fd883f3971a64de418759cf3f1e420371de62be634` | n/a | 24 OCPUs, 48 threads |
| 7 | v6.2 recomputation of segment 2 | 4,000,319,000,000,000,000 → 4.001·10^18 − 2 | the run 5 binary | `bfea4378691c270f0d12e1a4d65b1f1d059c74ce370ee0e7da723bc2b017a370` | `cd2d1cd0…e683` | 24 OCPUs, 48 threads |

PGO profiles were trained on one thread over [4·10^18, 4·10^18 + 2·10^8 − 2] with
block bits 24, anchor limit 4·10^7 and ring 512.

## Notes on the records

- **Run 2.** `provenance.txt` was written at build time and is covered by the
  campaign's own `SHA256SUMS`. The record repeats the segment 2 line, the second
  time with an empty `start=`. The segment 2 start is 4,000,319,000,000,000,000, as
  shown in the first line and in the segment 2 log and checkpoint.
- **Runs 3 and 4.** These records were reconstructed on 2026-10-01 from the build
  output printed on the VM at build time. The original files were lost with the VM.
  The source and binary hashes are the ones printed at build time. For run 3 the
  compiler version was not recorded. Run 4 was paused once for a benchmark and
  resumed with the same binary (`--resume`).
- **Run 5.** The campaign source was rebuilt on 2026-10-01 with the same settings,
  and the rebuild produced a binary byte-identical to the campaign binary
  (`bfea4378…a370`). The rebuild also passed the bit-exact SIMD check, the forced-miss
  comparison between the PGO and plain builds, and a row-by-row certificate dump at
  campaign settings (see `provenance/5_extension_v62_build_reproduction.txt`).
- **Run 6.** This run uses the v5 campaign source, built with a different compiler
  (GCC 13.4.0), so its binary differs in bytes from the run 2 binary. The binary is
  archived with the data.
- **Run 7.** This run used the run 5 campaign binary (`bfea4378…a370`), whose
  reproduction from the archived v6.2 source is described under run 5.
- **Campaign 1.** Only checkpoints survive: 13 files, with no log, build record or
  hardware record. The headers show three invocations that join end to end, each
  starting on the last number verified by the one before:
  [4·10^18, 4,000,037,279,999,999,998], then up to 4,000,121,119,999,999,996, then
  up to 4.001·10^18 (resumed several times, with a checkpoint saved at each stop).
  Every header has `total_verified` equal to the number of even integers in its range,
  and zero misses. The archived source `src/archive/goldbach_verifier_ckpt_be0158aa.cpp`
  writes checkpoints in this format and samples QHOT and SIEVE rows the way
  checkpoints 1–12 do. The final checkpoint is different: its 596 rows are all
  SIEVE, at spaced N, and 18 of them are not minimal partitions (their p is still
  prime and p + q = N). be0158aa cannot write rows like that. So the last stretch,
  from 4,000,914,939,999,999,996 to 4.001·10^18 (about 4.3·10^13 numbers), was
  probably run by a later revision whose source is not archived. Campaign 1 is not
  needed for the result: runs 2, 3 and 7 cover the same range with full records.
