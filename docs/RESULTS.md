# Results

Every even integer in [4·10^18, 4.003·10^18] was verified to be a sum of two primes
(1.5·10^15 numbers, zero misses). Each part of the range was verified by two
independent runs, both with full provenance (source, binary, log, checkpoint). The
first 5·10^14 numbers were also verified a third time, by campaign 1. Run numbers
match the folders of the Zenodo data record and [`PROVENANCE.md`](PROVENANCE.md).

## Coverage

| Range | Primary run | Independent rerun | Also verified by | COLD witnesses compared | Identical |
| --- | --- | --- | --- | --- | --- |
| [4.000, 4.001)·10^18 | 2: v5 | 3: v6 (segment 1) and 7: v6.2 (segment 2) | 1: campaign 1 | 100,000 | all |
| [4.001, 4.002)·10^18 | 4: v6, K=40 | 6: v5 | | 100,000 | all |
| [4.002, 4.003)·10^18 | 5: v6.2, K=56 | 6: v5 | | 100,000 | all |

COLD witnesses are minimal partitions, so two correct runs with the same block size
must produce the same COLD rows. QHOT rows may legitimately differ. Here 57–58% of
QHOT rows coincide, because the two programs use different rings.

## Runs

| | 2: v5 seg. 1 | 2: v5 seg. 2 | 3: v6 rerun seg. 1 | 7: v6.2 rerun seg. 2 | 4: v6 ext. | 5: v6.2 ext. | 6: v5 rerun ext. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| First N | 4,000,000,000,000,000,000 | 4,000,319,000,000,000,000 | 4,000,000,000,000,000,000 | 4,000,319,000,000,000,000 | 4,001,000,000,000,000,000 | 4,002,000,000,000,000,000 | 4,001,000,000,000,000,000 |
| Last N | 4,000,318,999,999,999,998 | 4,000,999,999,999,999,998 | 4,000,318,999,999,999,998 | 4,000,999,999,999,999,998 | 4,001,999,999,999,999,998 | 4,002,999,999,999,999,998 | 4,002,999,999,999,999,998 |
| Even integers | 1.595·10^14 | 3.405·10^14 | 1.595·10^14 | 3.405·10^14 | 5·10^14 | 5·10^14 | 10^15 |
| Anchor limit / ring | 10^5 / 64 | 10^6 / 64 | 4·10^7 / 512 | 4·10^7 / 512 | 4·10^7 / 512 | 4·10^7 / 512 | 10^6 / 64 |
| Misses | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Ring hit rate | 99.916694% | 99.983884% | 99.999132% | 99.999132% | 99.999132% | 99.999132% | 99.983884% |
| Miller–Rabin calls per N | 0.00644 | 0.00124 | 0.000070 | 0.000070 | 0.000070 | 0.000070 | 0.00124 |
| Wall time | 16.5 h | 17.2 h | 2.39 h | 2.28 h | 7.25 h | 3.36 h | 25.5 h |
| Throughput (M/s) | 2,689 | 5,494 | 18,561 | 41,494 | 19,145 | 41,373 | 10,874 |
| VM | 12 OCPUs | 12 OCPUs | 24 OCPUs | 24 OCPUs | 24 OCPUs | 24 OCPUs | 24 OCPUs |
| Witness rows | 63,800 | 136,200 | 63,800 | 136,200 | 200,000 | 200,000 | 400,000 |

All block bits are 24. Hardware for runs 2–7: OCI VM.Standard.E5.Flex, AMD EPYC
9J14, two threads per OCPU. Run 4 was paused once for a benchmark and resumed with
the same binary, and its wall time is the total of both parts. Segment 1 of run 2
was stopped at 4,000,318,999,999,999,998, and segment 2 started at the next even
number with a larger anchor limit.

## Campaign 1 (run 1)

Campaign 1 used an earlier verifier with per-number Q-Hot and a segmented-sieve
prefilter ahead of Miller–Rabin. It ran as three invocations that join end to end:

| Invocation | First N | Last N | Even integers | Misses | Ring hit rate | Wall time | M/s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4,000,000,000,000,000,000 | 4,000,037,279,999,999,998 | 1.864·10^13 | 0 | 99.986374% | 5.5 h | 939 |
| 2 | 4,000,037,279,999,999,998 | 4,000,121,119,999,999,996 | 4.192·10^13 | 0 | 99.986374% | 12.4 h | 936 |
| 3 | 4,000,121,119,999,999,996 | 4,001,000,000,000,000,000 | 4.394·10^14 | 0 | 99.961466% | 130.3 h | 937 |

Each invocation starts on the last number of the one before, so the campaign covers
every even integer in [4·10^18, 4.001·10^18], both ends included. In every
checkpoint header, `total_verified` equals the number of even integers in its range.
It ran on an OCI E5.Flex VM with 12 OCPUs (AMD EPYC 9J14), according to the author's records. Only the checkpoints survive, with no logs or build records. See
[`PROVENANCE.md`](PROVENANCE.md) for what this means for the final stretch of
invocation 3.

## Witness checks

Every witness row of runs 2–7 (1,200,000 rows) was re-checked with
`scripts/verify_witnesses.py --minimal` (SymPy 1.14). The checks were: N even,
p + q = N, p and q prime, p within the anchor limit, N inside the verified range, and
each COLD witness minimal. All rows passed. The largest COLD p is 2,237, well inside
the script's minimality search bound of 20,000, so the minimality check is complete.

The 4,540,596 distinct witness rows of campaign 1 (QHOT and SIEVE, from all 13
checkpoints) were checked the same way. All are valid partitions: p + q = N, with p
and q prime. Two SIEVE rows fall on the same N as a COLD row of a later run, and they
agree. 18 SIEVE rows, all in the final checkpoint, are not minimal. These rows still
certify their N, but they show that the final checkpoint was not written by the
archived source.

## Performance on one 12-OCPU EPYC 9J14 VM (24 threads, same work)

| Build | M evens/s | vs v4 |
| --- | --- | --- |
| v4 (per-N QHot) | 915 | 1.0× |
| v5 (batched), anchor 10^5 | 2,689 | 2.9× |
| v5, anchor 10^6 | 5,494 | 6.0× |
| v6 sorted ring + PGO (bb 24, 4·10^7, ring 512) | about 7,100 | 7.8× |
| v6 + 32-probe prefix | 8,280 | 9.0× |
| v6 + 32-probe prefix + AVX-512 | 9,331 | 10.2× |
| v6 + 40-probe AVX-512 prefix | 9,618 | 10.5× |
| v6.1: + fully-covered-batch shortcut | 14,825 | 16.2× |
| **v6.2: + carried words in the AVX-512 prefix** | **18,641** | **20.4×** |

24 OCPUs (48 threads): 19,214 M/s for v6 (K=40, AVX-512) and 29,580 M/s for
v6.1 (about 2.1 cycles per even number per core at the reported 2.6 GHz).
v6.2 on 12 OCPUs: 18,641 M/s; on 24 OCPUs (48 threads): 37,094 M/s, with v6 at
19,138 and v6.1 at 29,525 in the same session (K=40).
With carried words the best prefix length rose from 40 to 56: v6.2 with K=56 runs at
41,526 M/s on 24 OCPUs (K=40: 37,185; K=48: 40,295; K=64: 38,894), about 1.5 cycles
per even number per core. Thread pinning (OMP_PROC_BIND) made no difference.

## Cross-platform (v6 final settings, PGO trained per machine)

| Machine | ISA | Threads | Result | Hit rate |
| --- | --- | --- | --- | --- |
| EPYC 9J14, 24 OCPUs | AVX-512 | 48 | 18,547 M/s | 99.999133% |
| EPYC 9J14, 12 OCPUs | AVX-512 | 24 | 9,331 M/s | 99.999133% |
| Apple M4 (MacBook Air) | NEON | 10 | about 6,200 M/s burst, about 3,800 sustained | 99.999132% |
| Intel Alder Lake laptop | AVX2 | 12 | thermally limited | 99.999132% |

Hit rates agree to six decimals across architectures (differences in the last
digit reflect different measured ranges).
