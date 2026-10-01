#!/usr/bin/env bash
# Build and check a native v6.2 binary on macOS (Apple Silicon, Homebrew GCC, python3 + SymPy).
#   - uses the NEON prefix with K=56 (measured best on Apple M4; carried words off by default)
#   - trains a PGO profile at exactly the campaign settings (PGO gave +78-86% on the M4)
#   - runs the bit-exact prefix check, the forced-miss check and a certificate dump
#   - writes build/provenance_v6.txt (settings, source, binary and profile SHA-256)
# Usage (from the repository root, in bash): bash scripts/build_v6_mac.sh
# Optional overrides: GB_CXX=g++-15 (compiler), KFIX=n (prefix length, a multiple of 4)
set -e
SRC=${1:-src/goldbach_v6_2.cpp}
BB=24; P=40000000; RING=512; TR=200000000; S=4000000000000000000

# --- tools ---------------------------------------------------------------
CXX=${GB_CXX:-}
if [ -z "$CXX" ]; then
  for v in 17 16 15 14 13 12 11; do
    if command -v g++-$v > /dev/null 2>&1; then CXX=g++-$v; break; fi
  done
fi
[ -n "$CXX" ] || { echo "No Homebrew GCC found (g++-NN). Install it with: brew install gcc"; exit 1; }
"$CXX" --version 2>/dev/null | head -1 | grep -qi clang && { echo "$CXX is Apple clang, not GCC. Install GCC with: brew install gcc"; exit 1; }
python3 -c "import sympy" 2>/dev/null || { echo "SymPy missing. Install it with: pip3 install --user sympy"; exit 1; }
T=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
if command -v shasum > /dev/null 2>&1; then SHA="shasum -a 256"; else SHA="sha256sum"; fi
CPU=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || grep -m1 "model name" /proc/cpuinfo | cut -d: -f2)

F=${GB_FLAGS:-"-O3 -mcpu=native -fopenmp -std=c++17"}
ISA=${GB_ISA--DSIMD_NEON}
K=${KFIX:-56}
[ $((K % 4)) -eq 0 ] || { echo "KFIX must be a multiple of 4 for NEON"; exit 1; }
KF="-DKFIX=$K"

mkdir -p build && cp "$SRC" build/v6.cpp && cp scripts/verify_dump.py build/ && cd build
echo "== compiler: $CXX | prefix: ${ISA:-scalar}, K=$K | threads: $T"
"$CXX" $F $ISA $KF v6.cpp -o goldbach_v6_plain
"$CXX" $F $KF -DDUMP_ALL v6.cpp -o goldbach_v6_dump
rm -f gs.o gs.gcda
"$CXX" $F $ISA $KF -fprofile-generate -c v6.cpp -o gs.o && "$CXX" $F -fprofile-generate gs.o -o gs_gen
mkdir -p train && (cd train && rm -f *.csv; ../gs_gen $S $((S+TR-2)) $BB 1 $P 1000 $RING > /dev/null)
[ -s gs.gcda ] || { echo "PGO profile missing (gs.gcda) - STOP"; exit 1; }
"$CXX" $F $ISA $KF -fprofile-use -Wno-missing-profile -c v6.cpp -o gs.o && "$CXX" $F gs.o -o goldbach_v6

if [ -n "$ISA" ]; then
  echo "== bit-exact check of the vector prefix (K=$K)"
  "$CXX" $F $ISA $KF -DSIMD_CHECK v6.cpp -o goldbach_v6_check
  mkdir -p chk_simd && (cd chk_simd && rm -f *.csv
    out=$(../goldbach_v6_check $S $((S+2*16777216-2)) $BB $T $P 1000 $RING 2>&1 | grep -E "SIMD checks|MISMATCH" || true)
    echo "$out"; echo "$out" | grep -q "SIMD checks passed" || { echo "SIMD CHECK FAILED - STOP"; exit 1; })
fi

echo "== forced misses: PGO and plain builds must agree"
for b in goldbach_v6 goldbach_v6_plain; do
  mkdir -p chk_$b && (cd chk_$b && rm -f *.csv; ../$b $S $((S+2*4194304-2)) 16 $T 1200 1000 $RING > out.txt 2> err.txt || true)
done
grep -o "MISS.*" chk_goldbach_v6/err.txt > chk_m1.txt || true
grep -o "MISS.*" chk_goldbach_v6_plain/err.txt > chk_m2.txt || true
if [ -s chk_m1.txt ] && cmp -s chk_m1.txt chk_m2.txt; then
  echo "MISS LISTS IDENTICAL ($(wc -l < chk_m1.txt | tr -d ' ') misses)"
else echo "MISS LISTS DIFFER OR EMPTY - STOP"; exit 1; fi

echo "== certificate dump at campaign settings (2^20 evens)"
mkdir -p chk_dump && cd chk_dump && rm -f *.csv
../goldbach_v6_dump $S $((S+2*1048576-2)) $BB $T $P 1000 $RING > /dev/null
python3 ../verify_dump.py goldbach_v62_dump.csv $S $((S+2*1048576-2)) $P
rm -f goldbach_v62_dump.csv; cd ..

{ date -u; echo "CPU: $CPU"; echo "threads=$T"; "$CXX" --version | head -1
  echo "prefix: ${ISA:-scalar}; KFIX=$K; campaign settings: bb=$BB P=$P ring=$RING"
  echo "compile flags: $F $ISA $KF (+ PGO)"
  echo "PGO training: 1 thread, [$S, $((S+TR-2))], bb=$BB P=$P ring=$RING"
  $SHA v6.cpp goldbach_v6 gs.gcda; } | tee provenance_v6.txt
echo "== done: build/goldbach_v6"
echo "   run: ./build/goldbach_v6 START END $BB $T $P 250000 $RING"
