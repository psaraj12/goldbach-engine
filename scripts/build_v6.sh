#!/usr/bin/env bash
# Build and check a native v6 binary on this machine (Linux, GCC >= 11, python3-sympy).
#   - picks the vector prefix: AVX-512 > AVX2 > NEON (ARM64) > scalar
#   - picks the prefix length K for that prefix (override with KFIX=n)
#       AVX-512 and NEON: 56 (EPYC 9J14, Apple M4); AVX2 and scalar: 32
#   - trains a PGO profile at exactly the campaign settings
#   - runs the bit-exact prefix check, the forced-miss check and a certificate dump
#   - writes build/provenance_v6.txt (settings, source, binary and profile SHA-256)
# Usage (from the repository root): bash scripts/build_v6.sh [path/to/source.cpp]
#   (default src/goldbach_v6_3.cpp)
set -e
SRC=${1:-src/goldbach_v6_3.cpp}
# Settings follow the source version: v6.3 uses anchor limit 2e7 and, with PGO, the carried
# index (-DCARRYIDX); earlier versions use 4e7 and no extra flag, so the archived v6.2
# campaign source still rebuilds exactly as it did. Overrides: GB_P=n, GB_EXTRA="flags".
if grep -q 'PROG_VERSION = "goldbach_v6_3' "$SRC"; then PDEF=20000000; XDEF="-DCARRYIDX"
else PDEF=40000000; XDEF=""; fi
P=${GB_P:-$PDEF}; X=${GB_EXTRA-$XDEF}
BB=24; RING=512; TR=200000000; S=4000000000000000000; T=$(nproc)
if [ "$(uname -m)" = aarch64 ]; then
  F="-O3 -mcpu=native -fopenmp -std=c++17"; ISA=-DSIMD_NEON; KDEF=56
else
  F="-O3 -march=native -fopenmp -std=c++17"
  if   grep -qw avx512f /proc/cpuinfo; then ISA=-DSIMD_AVX512; KDEF=56
  elif grep -qw avx2    /proc/cpuinfo; then ISA=-DSIMD_AVX2;   KDEF=32
  else ISA="";                                   KDEF=32; fi
fi
K=${KFIX:-$KDEF}
case "$ISA" in
  -DSIMD_AVX512) [ $((K % 8)) -eq 0 ] || { echo "KFIX must be a multiple of 8 for AVX-512"; exit 1; } ;;
  -DSIMD_AVX2|-DSIMD_NEON) [ $((K % 4)) -eq 0 ] || { echo "KFIX must be a multiple of 4 for $ISA"; exit 1; } ;;
esac
KF="-DKFIX=$K"
mkdir -p build && cp "$SRC" build/v6.cpp && cp scripts/verify_dump.py build/ && cd build
echo "== vector prefix: ${ISA:-scalar}, K=$K, anchor limit $P${X:+, $X}"
g++ $F $ISA $KF $X v6.cpp -o goldbach_v6_plain
g++ $F $KF $X -DDUMP_ALL v6.cpp -o goldbach_v6_dump
rm -f gs.o gs.gcda
g++ $F $ISA $KF $X -fprofile-generate -c v6.cpp -o gs.o && g++ $F -fprofile-generate gs.o -o gs_gen
mkdir -p train && (cd train && rm -f *.csv; ../gs_gen $S $((S+TR-2)) $BB 1 $P 1000 $RING > /dev/null)
[ -s gs.gcda ] || { echo "PGO profile missing - STOP"; exit 1; }
g++ $F $ISA $KF $X -fprofile-use -Wno-missing-profile -c v6.cpp -o gs.o && g++ $F gs.o -o goldbach_v6
if [ -n "$ISA" ]; then
  echo "== bit-exact check of the vector prefix (K=$K)"
  g++ $F $ISA $KF $X -DSIMD_CHECK v6.cpp -o goldbach_v6_check
  mkdir -p chk_simd && (cd chk_simd && rm -f *.csv
    out=$(../goldbach_v6_check $S $((S+2*33554432-2)) $BB $T $P 1000 $RING 2>&1 | grep -E "SIMD checks|MISMATCH")
    echo "$out"; echo "$out" | grep -q "SIMD checks passed" || { echo "SIMD CHECK FAILED - STOP"; exit 1; })
fi
echo "== forced misses: PGO and plain builds must agree"
for b in goldbach_v6 goldbach_v6_plain; do
  mkdir -p chk_$b && (cd chk_$b && rm -f *.csv; ../$b $S $((S+2*4194304-2)) 16 $T 1200 1000 $RING > out.txt 2> err.txt || true)
done
diff <(grep -o "MISS.*" chk_goldbach_v6/err.txt) <(grep -o "MISS.*" chk_goldbach_v6_plain/err.txt) > /dev/null \
  && echo "MISS LISTS IDENTICAL ($(grep -c MISS chk_goldbach_v6/err.txt) misses)" || { echo "MISS LISTS DIFFER - STOP"; exit 1; }
echo "== certificate dump at campaign settings (2^20 evens)"
mkdir -p chk_dump && cd chk_dump && rm -f *.csv
../goldbach_v6_dump $S $((S+2*1048576-2)) $BB $T $P 1000 $RING > /dev/null
D=$(ls goldbach_v6*_dump.csv); python3 ../verify_dump.py $D $S $((S+2*1048576-2)) $P
rm -f $D; cd ..
{ date -u; lscpu 2>/dev/null | grep "Model name"; echo "threads=$T"; g++ --version | head -1
  echo "prefix: ${ISA:-scalar}; KFIX=$K; campaign settings: bb=$BB P=$P ring=$RING"
  echo "compile flags: $F $ISA $KF${X:+ $X} (+ PGO)"
  echo "PGO training: 1 thread, [$S, $((S+TR-2))], bb=$BB P=$P ring=$RING"
  sha256sum v6.cpp goldbach_v6 gs.gcda; } | tee provenance_v6.txt
echo "== done: build/goldbach_v6"
echo "   run: ./build/goldbach_v6 START END $BB \$(nproc) $P 250000 $RING"
