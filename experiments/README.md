# Experiments

Variants measured during development and kept for the record. None is used in
production; results are single-thread unless stated (4·10^18, same work).

| File | Idea | Outcome |
| --- | --- | --- |
| goldbach_v5_b32.cpp | 32-number batches | −8% (more probes per N) |
| goldbach_v5_b128.cpp | 128-number batches | −3% single thread, −7% on EPYC |
| goldbach_v5_desc.cpp | cold path scans anchors in descending order | about 70× slower (ring entries expire immediately) |
| goldbach_v5_tc.cpp | anchors restricted to twin/cousin primes | about 3× slower (fewer partners, clustered coverage) |
| goldbach_v5_mod6.cpp | ring split by q mod 6, walked alternately | no measurable gain (ring already mixes classes) |
| goldbach_v5_sorted.cpp | ring sorted by q, early exit | +5–15% (adopted in v6) |
| goldbach_v5_kfix.cpp | fixed unrolled prefix of K probes | +16% on EPYC at K=32 (adopted in v6); with AVX-512, K=40 is +3.6% more |
| (in src/goldbach_v6_1.cpp) | skip the per-number loop when the ring covers the whole batch | +54% on EPYC, 24 and 48 threads, identical output (adopted in v6.1) |

Four independent accumulators in the prefix (`-DKACC4` in `src/goldbach_v6_simd.cpp`)
were −3% on x86 and −12% on Apple M4, and are kept only as a documented option.
| goldbach_v6_fs.cpp | fully-covered-batch shortcut as a flag (-DFULLSKIP) | +54% on EPYC, identical output (adopted in v6.1) |
| goldbach_v6_1_carry.cpp (-DCARRY) | carry each prefix entry's high word to the next batch | +25.7% on EPYC 12 OCPUs, identical output (adopted in v6.2) |
| goldbach_v6_1_carry.cpp (-DVTAIL) | vectorize the walk after the prefix, 8 entries per step | −3.4% on EPYC (extra probes and gathers outweigh the saving) |
| goldbach_v6_2_cw.cpp (-DCLASSWALK) | after the prefix, skip ring entries whose class mod 3 cannot fill the remaining gaps | probes −6% but 12–20% slower (per-entry mod 3 and unpredictable branch) |
| (in src/goldbach_v6_2.cpp, -DCARRY) | carried words for the scalar prefix | 5–16% slower on x86; off by default for scalar |
| (in src/goldbach_v6_2.cpp, -DKFIX) | prefix length re-tuned after carried words: K = 40/48/56/64 | 37,185 / 40,295 / 41,526 / 38,894 M/s on 24 OCPUs; K=56 adopted |
| (in src/goldbach_v6_2.cpp, -DCKPT_EVENS) | checkpoint chunk 10^11 / 2^37 / 2^40 evens | within 0.2%; default 10^11 kept (chunk size also fixes block layout, so keep it for comparable runs) |
| (environment) | OMP_PROC_BIND / OMP_PLACES thread placement | no effect (within ±0.2%) |
