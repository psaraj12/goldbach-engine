"""Independent audit of a Conjecture-A (twin/cousin-starter) campaign.
Usage: python3 audit_tc.py <results_dir>
Expects goldbach_tc_checkpoint.csv, goldbach_tc_misses.csv, goldbach_tc_minp.csv.
Rechecks every sampled witness and every logged minimal witness with sympy,
spot-checks minimality and the three per-gear minima, and verifies the mod-3
gearbox structure (Proposition 3 / Corollary 4 of the paper)."""
import os, random, statistics, sys
from sympy import isprime

D = sys.argv[1] if len(sys.argv) > 1 else "."
def starter(p):  # twin or cousin starter
    return isprime(p) and (isprime(p + 2) or isprime(p + 4))

hdr, samples = {}, []
for line in open(os.path.join(D, "goldbach_tc_checkpoint.csv")):
    if line.startswith("# ") and "=" in line:
        k, v = line[2:].strip().split("=", 1); hdr[k] = v
    elif line[0].isdigit():
        N, p, q, src = line.strip().split(","); samples.append((int(N), int(p), int(q), src))
lo, hi = int(hdr["original_start"]), int(hdr["original_end"])
expect = (hi - lo) // 2 + 1
P = int(hdr["cfg_p_anchor_limit"])
print("CHECKPOINT")
print(f"  version={hdr['cfg_version']}  range=[{lo},{hi}]  last_verified={hdr['last_verified']}  P={P}")
print(f"  total_verified={hdr['total_verified']}  expected evens={expect}  "
      f"match={int(hdr['total_verified']) == expect == int(hdr['total_total'])}")
print(f"  total_misses={hdr['total_misses']} unresolved={hdr['unresolved']}")
print(f"  ring hit rate={100 * int(hdr['total_qhot']) / expect:.6f}%  "
      f"mr_calls==cold_survivors: {hdr['mr_calls'] == hdr['cold_survivors']}")
bad = [s for s in samples if not (s[1] + s[2] == s[0] and starter(s[1]) and isprime(s[2])
                                  and lo <= s[0] <= hi and s[1] <= P)]
print(f"  samples={len(samples)}  all witnesses valid (p twin/cousin starter, q prime, p+q=N): {not bad}  bad={bad[:3]}")

mrows = [l for l in open(os.path.join(D, "goldbach_tc_misses.csv")) if l[0].isdigit()]
print(f"MISSES file: {len(mrows)} rows")

rows = []
for line in open(os.path.join(D, "goldbach_tc_minp.csv")):
    if line[0].isdigit():
        N, p, q, t, pt, pc = line.strip().split(",")
        rows.append((int(N), int(p), int(q), t, int(pt), int(pc)))
print(f"MINP rows={len(rows)}  N range [{min(r[0] for r in rows):,} .. {max(r[0] for r in rows):,}]")
badm = sum(1 for N, p, q, *_ in rows if p + q != N or not (starter(p) and isprime(q)))
print(f"  rows with invalid witness: {badm}")
# type column and per-gear minima consistency
def tcode(p):
    t = (1 if isprime(p + 2) else 0) | (2 if isprime(p + 4) else 0)
    return {1: "T", 2: "C", 3: "B"}[t]
badt = sum(1 for r in rows[:50000] if tcode(r[1]) != r[3])
print(f"  type column matches (first 50,000 rows): mismatches = {badt}")
badg = sum(1 for r in rows if r[1] != min(x for x in (r[4], r[5]) if x))
print(f"  p == min(nonzero p_twin, p_cousin) for every row: {badg == 0}")
# gearbox: N=1 mod 3 never needs +4, N=0 never needs +4; N=2 needs +4 iff N-3 composite
g_ok = all(r[4] != 0 for r in rows if r[0] % 3 in (0, 1))
g2 = [r for r in rows if r[0] % 3 == 2]
g2_ok = all((r[4] == 0) == (not isprime(r[0] - 3)) for r in g2)
print(f"  classes N=0,1 (mod 3) always have a twin partner: {g_ok}")
print(f"  class N=2: twin partner exists exactly when N-3 is prime (Prop. 3): {g2_ok}")
# minimality spot-check
random.seed(1)
st = [p for p in range(3, 20000) if starter(p)]
tw = [p for p in st if isprime(p + 2)]; co = [p for p in st if isprime(p + 4)]
nonmin = 0
for N, p, q, t, pt, pc in random.sample(rows, 2000):
    if next(s for s in st if isprime(N - s)) != p: nonmin += 1
    if pt and next(s for s in tw if isprime(N - s)) != pt: nonmin += 1
    if pc and next(s for s in co if isprime(N - s)) != pc: nonmin += 1
print(f"  minimality spot-check (2000 rows, overall/twin/cousin minima): wrong = {nonmin}")
ps = [r[1] for r in rows]; mx = max(rows, key=lambda r: r[1])
print(f"  least starter: mean {statistics.fmean(ps):.1f}  max {mx[1]:,} at N = {mx[0]:,}")
