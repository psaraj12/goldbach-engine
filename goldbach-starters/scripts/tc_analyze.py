#!/usr/bin/env python3
# tc_analyze.py -- offline analysis of goldbach_tc_minp.csv (minp_format=2).
#
# Reads the [TC-INST] minimal-partner log written by goldbach_v6_4_tc and
# reports the research quantities of the twin/cousin propagation framework:
#   * gear usage by N mod 3 (the mod-3 complementarity mechanism)
#   * the needs-+4 fraction (p_twin == 0), tested against the 1/3 limit law
#   * growth of the minimal starter (fuel gauge): max / mean / median / p99
#   * log2 histogram of the minimal starter
#
# Usage: python3 tc_analyze.py [goldbach_tc_minp.csv]
import sys, statistics

path = sys.argv[1] if len(sys.argv) > 1 else "goldbach_tc_minp.csv"
rows = []
legacy = 0
with open(path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("N,"):
            continue
        parts = line.split(",")
        if len(parts) == 6:
            N, p, q, ty, ptw, pco = parts
            rows.append((int(N), int(p), int(q), ty, int(ptw), int(pco)))
        else:
            legacy += 1  # old 3-column rows carry no gear data

if not rows:
    sys.exit(f"no format-2 rows found in {path}" +
             (f" ({legacy} legacy rows skipped)" if legacy else ""))
if legacy:
    print(f"NOTE: {legacy} legacy 3-column row(s) skipped (no gear data)\n")

rows.sort()
Ns     = [r[0] for r in rows]
pmins  = [r[1] for r in rows]
n      = len(rows)

print(f"rows: {n:,}   N range: [{Ns[0]:,} .. {Ns[-1]:,}]")

# ---- gear usage by N mod 3 ----------------------------------------------
print("\n== Mod-3 gearbox (per residue class of N) ==")
print(f"{'N mod 3':>8} {'count':>9} {'needs +4':>9} {'%':>7}   twin-min  cousin-min  both-min")
tot_nc = 0
for m in (0, 1, 2):
    cls = [r for r in rows if r[0] % 3 == m]
    if not cls:
        continue
    nc  = sum(1 for r in cls if r[4] == 0)          # no twin partner <= limit
    tot_nc += nc
    tT  = sum(1 for r in cls if r[3] == 'T')
    tC  = sum(1 for r in cls if r[3] == 'C')
    tB  = sum(1 for r in cls if r[3] == 'B')
    print(f"{m:>8} {len(cls):>9,} {nc:>9,} {100*nc/len(cls):>6.2f}%"
          f"   {tT:>8,} {tC:>10,} {tB:>8,}")
print(f"\nneeds-+4 fraction overall: {100*tot_nc/n:.3f}%   (limit law predicts -> 33.333% as N -> inf)")

# ---- fuel gauge ----------------------------------------------------------
print("\n== Minimal starter (fuel gauge) ==")
mx = max(rows, key=lambda r: r[1])
srt = sorted(pmins)
p99 = srt[int(0.99 * (n - 1))]
print(f"max p_min   : {mx[1]:,}   at N = {mx[0]:,}")
print(f"mean / median / p99 : {statistics.fmean(pmins):,.1f} / {statistics.median(pmins):,.0f} / {p99:,}")
tw = [r[4] for r in rows if r[4]]
co = [r[5] for r in rows if r[5]]
if tw: print(f"minimal TWIN partner   (where it exists): max {max(tw):,}, median {statistics.median(tw):,.0f}  [{len(tw):,} rows]")
if co: print(f"minimal COUSIN partner (where it exists): max {max(co):,}, median {statistics.median(co):,.0f}  [{len(co):,} rows]")

# ---- per-stretch crawl curve ---------------------------------------------
# Windows fire once per 10^10 stretch, so grouping rows by N // 10^10 gives
# one point per window: the needs-+4 fraction and max minimal starter along
# the range. Printed only when the data spans more than one stretch; also
# written to tc_crawl.csv for plotting.
STRETCH = 10**10
stretches = {}
for r in rows:
    s = r[0] // STRETCH
    c, nc, mp = stretches.get(s, (0, 0, 0))
    stretches[s] = (c + 1, nc + (1 if r[4] == 0 else 0), max(mp, r[1]))
if len(stretches) > 1:
    print("\n== Crawl curve (one point per 10^10 stretch) ==")
    with open("tc_crawl.csv", "w") as out:
        out.write("stretch_start,rows,needs4_pct,max_p_min\n")
        for s in sorted(stretches):
            c, nc, mp = stretches[s]
            out.write(f"{s*STRETCH},{c},{100*nc/c:.4f},{mp}\n")
    ss = sorted(stretches)
    show = ss if len(ss) <= 12 else ss[:4] + ss[len(ss)//2-2:len(ss)//2+2] + ss[-4:]
    print(f"{'stretch start':>16} {'rows':>7} {'needs +4':>9} {'max p_min':>10}")
    prev = None
    for s in show:
        if prev is not None and s != prev + 1:
            print(f"{'...':>16}")
        c, nc, mp = stretches[s]
        print(f"{s*STRETCH:>16,} {c:>7,} {100*nc/c:>8.2f}% {mp:>10,}")
        prev = s
    first, last = stretches[ss[0]], stretches[ss[-1]]
    print(f"full curve ({len(ss)} points) -> tc_crawl.csv   "
          f"[first stretch {100*first[1]/first[0]:.2f}% .. last {100*last[1]/last[0]:.2f}%]")

# ---- histogram -----------------------------------------------------------
print("\n== p_min histogram (log2 buckets) ==")
hist = {}
for p in pmins:
    b = p.bit_length()          # bucket: 2^(b-1) <= p < 2^b
    hist[b] = hist.get(b, 0) + 1
peak = max(hist.values())
for b in sorted(hist):
    c = hist[b]
    bar = "#" * max(1, round(40 * c / peak))
    print(f"[{1 << (b-1):>6,} .. {(1 << b) - 1:>6,}] {c:>8,}  {bar}")
