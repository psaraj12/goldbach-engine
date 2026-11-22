#!/usr/bin/env bash
# Build and check a native v6 binary on this machine (Linux, GCC >= 11, python3-sympy).
#   - picks the vector prefix: AVX-512 > AVX2 > NEON (ARM64) > scalar
#   - picks the prefix length K for that prefix (override with KFIX=n)
#       AVX-512 and NEON: 56 (EPYC 9J14, Apple M4); AVX2 and scalar: 32
#   - trains a PGO profile at exactly the campaign settings
#   - runs the bit-exact prefix check, the forced-miss check and a certificate dump
#   - writes <build dir>/provenance_v6.txt (settings, source, binary and profile SHA-256)
# Usage (from the repository root): bash scripts/build_v6.sh [path/to/source.cpp]
#   (default goldbach_v6_4.cpp)
#
# Anchor modes (selected by the source and GB_SEXY):
#   all odd primes        any source except goldbach_v6_4_tc* (the Goldbach campaign)
#   twin/cousin starters  goldbach_v6_4_tc*.cpp            (Conjecture A)
#   sexy starters         goldbach_v6_4_tc*.cpp, GB_SEXY=1  (Conjecture B, -DSEXY_ANCHORS)
#   e.g.  GB_SEXY=1 GB_BUILD_DIR=~/gb_build_sx bash scripts/build_v6.sh src/goldbach_v6_4_tc.cpp
#
# Overrides: GB_P=n (campaign anchor limit), GB_EXTRA="flags", GB_MISSP=n (forced-miss
#   anchor limit), GB_BUILD_DIR=dir (default build), GB_SEXY=1 (sexy-starter anchors).
set -e
SRC=${1:-goldbach_v6_4.cpp}
# Settings follow the source version: v6.3 uses anchor limit 2e7 and, with PGO, the carried
# index (-DCARRYIDX); earlier versions use 4e7 and no extra flag, so the archived v6.2
# campaign source still rebuilds exactly as it did.
# Starter-anchor builds (twin/cousin or sexy) keep the campaign limit 4e7 but need a
# HIGHER forced-miss limit, because the starter table is ~4x thinner than the prime
# table. Measured at 4e18 over 2^22 evens with limit 8000: 63 misses (twin/cousin,
# 345 starters) and 64 misses (sexy, 344 starters), comparable to the unrestricted
# build at limit 1200. Limit 1200 with starters forces ~150k misses and is very slow.
if grep -q 'PROG_VERSION = "goldbach_v6_4_tc' "$SRC"; then
  PDEF=40000000; XDEF=""; MISSPDEF=8000
  if [ "${GB_SEXY:-0}" = 1 ]; then
    grep -q 'SEXY_ANCHORS' "$SRC" || { echo "GB_SEXY=1 but $SRC has no SEXY_ANCHORS switch - STOP"; exit 1; }
    XDEF="-DSEXY_ANCHORS"; ANCHORS="sexy starters (p+6 prime)"
  else
    ANCHORS="twin/cousin starters (p+2 or p+4 prime)"
  fi
elif [ "${GB_SEXY:-0}" = 1 ]; then
  echo "GB_SEXY=1 needs a goldbach_v6_4_tc source - STOP"; exit 1
elif grep -q 'PROG_VERSION = "goldbach_v6_3' "$SRC"; then
  PDEF=20000000; XDEF="-DCARRYIDX"; MISSPDEF=1200; ANCHORS="all odd primes"
else
  PDEF=40000000; XDEF=""; MISSPDEF=1200; ANCHORS="all odd primes"
fi
P=${GB_P:-$PDEF}; X=${GB_EXTRA-$XDEF}; MISSP=${GB_MISSP:-$MISSPDEF}
# GB_EXTRA replaces the default extra flags; keep -DSEXY_ANCHORS if sexy mode was asked for.
if [ "${GB_SEXY:-0}" = 1 ] && [[ "$X" != *SEXY_ANCHORS* ]]; then X="$X -DSEXY_ANCHORS"; fi
BB=24; RING=512; TR=200000000; S=4000000000000000000; T=$(nproc)
# On WSL, a repository under /mnt/* pays 10-50x on the fsync-heavy check stages;
# set GB_BUILD_DIR=~/gb_build (native ext4) to build and test there.
BD=${GB_BUILD_DIR:-build}
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
SCRIPTS_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$BD" && cp "$SRC" "$BD"/v6.cpp && cp "$SCRIPTS_DIR"/verify_dump.py "$BD"/ && cd "$BD"
echo "== vector prefix: ${ISA:-scalar}, K=$K, anchor limit $P (anchors: $ANCHORS)${X:+, $X}"
g++ $F $ISA $KF $X v6.cpp -o goldbach_v6_plain
g++ $F $KF $X -DDUMP_ALL v6.cpp -o goldbach_v6_dump
rm -f gs.o gs.gcda
g++ $F $ISA $KF $X -fprofile-generate -c v6.cpp -o gs.o && g++ $F -fprofile-generate gs.o -o gs_gen
mkdir -p train && (cd train && rm -f *.csv; ../gs_gen $S $((S+TR-2)) $BB 1 $P 1000 $RING > /dev/null)
[ -s gs.gcda ] || { echo "PGO profile missing - STOP"; exit 1; }
g++ $F $ISA $KF $X -fprofile-use -Wno-missing-profile -c v6.cpp -o gs.o && g++ $F gs.o -o goldbach_v6
# Guard: the built binary must report the anchor mode that was asked for.
# (the usage text names the checkpoint file, which differs per anchor mode)
VER=$(./goldbach_v6 2>&1 | grep -o 'goldbach_[a-z0-9]*_checkpoint\.csv' | head -1 || true)
case "$ANCHORS" in
  sexy*)        [ "$VER" = goldbach_sx_checkpoint.csv ] || { echo "binary uses '$VER', expected goldbach_sx_checkpoint.csv - STOP"; exit 1; } ;;
  twin/cousin*) [ "$VER" = goldbach_tc_checkpoint.csv ] || { echo "binary uses '$VER', expected goldbach_tc_checkpoint.csv - STOP"; exit 1; } ;;
esac
if [ -n "$ISA" ]; then
  echo "== bit-exact check of the vector prefix (K=$K)"
  g++ $F $ISA $KF $X -DSIMD_CHECK v6.cpp -o goldbach_v6_check
  mkdir -p chk_simd && (cd chk_simd && rm -f *.csv
    out=$(../goldbach_v6_check $S $((S+2*33554432-2)) $BB $T $P 1000 $RING 2>&1 | grep -E "SIMD checks|MISMATCH")
    echo "$out"; echo "$out" | grep -q "SIMD checks passed" || { echo "SIMD CHECK FAILED - STOP"; exit 1; })
fi
echo "== forced misses (anchor limit $MISSP): PGO and plain builds must agree"
for b in goldbach_v6 goldbach_v6_plain; do
  mkdir -p chk_$b && (cd chk_$b && rm -f *.csv; ../$b $S $((S+2*4194304-2)) 16 $T $MISSP 1000 $RING > out.txt 2> err.txt || true)
done
diff <(grep -o "MISS.*" chk_goldbach_v6/err.txt) <(grep -o "MISS.*" chk_goldbach_v6_plain/err.txt) > /dev/null \
  && echo "MISS LISTS IDENTICAL ($(grep -c MISS chk_goldbach_v6/err.txt) misses)" || { echo "MISS LISTS DIFFER - STOP"; exit 1; }
echo "== certificate dump at campaign settings (2^20 evens)"
mkdir -p chk_dump && cd chk_dump && rm -f *.csv
../goldbach_v6_dump $S $((S+2*1048576-2)) $BB $T $P 1000 $RING > /dev/null
D=$(ls goldbach_v6*_dump.csv); python3 ../verify_dump.py $D $S $((S+2*1048576-2)) $P
rm -f $D; cd ..
{ date -u; lscpu 2>/dev/null | grep "Model name"; echo "threads=$T"; g++ --version | head -1
  echo "binary checkpoint file: $VER"
  echo "prefix: ${ISA:-scalar}; KFIX=$K; campaign settings: bb=$BB P=$P ring=$RING"
  echo "anchors: $ANCHORS; forced-miss limit: $MISSP"
  echo "compile flags: $F $ISA $KF${X:+ $X} (+ PGO)"
  echo "PGO training: 1 thread, [$S, $((S+TR-2))], bb=$BB P=$P ring=$RING"
  sha256sum v6.cpp goldbach_v6 gs.gcda; } | tee provenance_v6.txt
echo "== done: $BD/goldbach_v6"
case "$BD" in /*) RUNP="$BD" ;; *) RUNP="./$BD" ;; esac
echo "   run: $RUNP/goldbach_v6 START END $BB \$(nproc) $P 250000 $RING"
