# Results

## Campaign: [4·10^18, 4.001·10^18] (v5)

| | Segment 1 | Segment 2 | Campaign |
| --- | --- | --- | --- |
| Range | 4·10^18 → 4,000,318,999,999,999,998 | 4,000,319,000,000,000,000 → 4,000,999,999,999,999,998 | 4·10^18 → 4.001·10^18 − 2 |
| Even integers | 1.595·10^14 | 3.405·10^14 | 5.000·10^14 |
| Anchor limit / ring / block bits | 10^5 / 64 / 24 | 10^6 / 64 / 24 | |
| Exceptions | 0 | 0 | 0 |
| Ring hit rate | 99.9167% | 99.9839% | |
| Miller–Rabin calls per N | 0.00644 | 0.00124 | 0.00290 |
| Wall time | 16.5 h | 17.2 h | 33.7 h |
| Throughput (M/s) | 2,689 | 5,494 | 4,122 average |
| Witnesses | 63,800 | 136,200 | 200,000 |

Hardware: OCI VM.Standard.E5.Flex, 12 OCPUs (24 threads), AMD EPYC 9J14.
All 200,000 witnesses verified independently; all 100,000 COLD witnesses are
minimal partitions.

## Extension to 4.002·10^18 (v6)

[4.001·10^18, 4.002·10^18 − 2]: 5·10^14 even integers, zero exceptions, 7.25 h at
19,145 M/s on 24 OCPUs (v6, AVX-512, K=40). All 200,000 witnesses verified; all
100,000 COLD witnesses minimal. Contiguous with the v5 campaign.

## Independent recomputation (v6)

Both segments recomputed with v6 (AVX-512 prefix, anchor limit 4·10^7, ring 512)
on a 24-OCPU VM. Segment 1: zero exceptions; all 31,900 COLD witnesses identical to v5; all 63,800
v6 witnesses verified. Segment 2: *to be filled*.

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
