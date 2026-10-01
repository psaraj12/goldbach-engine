# Provenance

Each result is tied to an exact source file, binary and (for PGO builds) profile
by SHA-256. Fill in from the `provenance*.txt` files of each run before release.

| Run | Source file | Source SHA-256 | Binary SHA-256 | PGO profile SHA-256 | Hardware |
| --- | --- | --- | --- | --- | --- |
| v5 campaign, segment 1 | src/goldbach_v5_batch.cpp | *from provenance.txt* | *…* | n/a | OCI E5.Flex 12 OCPUs |
| v5 campaign, segment 2 | src/goldbach_v5_batch.cpp | *same as segment 1* | *…* | n/a | OCI E5.Flex 12 OCPUs |
| v6 rerun, segment 1 | src/goldbach_v6_simd.cpp | 6a7ace5fd6a37ba820a4b814067dc4f3582076f5d5f8b41d1ac7e46f5906f0af | 7d8f4ec3700effb6f6e8a2606028c85c82336077aee325566ba33c7d7c481efe | cc3890f2dfaf08b5d11fa7a985b6e10869e357f49844625d6c1ac5e543c592d7 | OCI E5.Flex 24 OCPUs |
| v6 extension 4.001→4.002·10^18 | src/goldbach_v6_simd.cpp (-DKFIX=40) | 6a7ace5f…f0af | 9d3cc280…301c | 4f7de34b…0687 | OCI E5.Flex 24 OCPUs |
| v6 rerun, segment 2 | src/goldbach_v6_simd.cpp | 6a7ace5f…f0af | *…* | *…* | OCI E5.Flex 24 OCPUs |

**Important:** the file `src/goldbach_v5_batch.cpp` in this repository must be
byte-identical to the campaign source (compare its SHA-256 with the campaign's
`provenance.txt`). If it differs, replace it with the archived campaign copy.
