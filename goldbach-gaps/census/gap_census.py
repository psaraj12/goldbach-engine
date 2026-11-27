"""Exhaustive g-starter Goldbach census.

For every gap g in GAPS and every even N with 6 <= N <= NMAX, find the least
g-starter p (p prime, p + g prime, p < N) with N - p prime, using ALL starters
below N (no anchor limit), so "exception" means no g-starter partition at all.

Output: one CSV row per gap.
"""
import sys, csv, time
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10**7
GAPS = range(6, 601, 6)
OUT  = sys.argv[2] if len(sys.argv) > 2 else "gap_census.csv"

GMAX = max(GAPS)
top = NMAX + GMAX + 2
s = np.ones(top + 1, dtype=bool); s[:2] = False
for i in range(2, int(top ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
primes = np.nonzero(s)[0]
odd_primes = primes[primes > 2]

N = np.arange(6, NMAX + 1, 2, dtype=np.int64)
twin_count = int(np.count_nonzero(s[odd_primes[odd_primes < NMAX] + 2]))

def hl_factor(g):
    """Hardy-Littlewood density of g-starters relative to twin starters."""
    f, m, l = 1.0, g, 3
    while m % 2 == 0: m //= 2
    while l * l <= m:
        if m % l == 0:
            f *= (l - 1) / (l - 2)
            while m % l == 0: m //= l
        l += 2
    if m > 1: f *= (m - 1) / (m - 2)
    return f

rows = []
for g in GAPS:
    t0 = time.time()
    cand = odd_primes[odd_primes < NMAX]
    st = cand[s[cand + g]]                     # g-starters below NMAX
    best = np.zeros(len(N), dtype=np.int64)
    idx = np.arange(len(N))                  # N not yet covered
    exc = []
    for p in st:                              # ascending: first hit = least starter
        if idx.size == 0: break
        done = N[idx] <= p                    # no starter < N left for these: exceptions
        if done.any():
            exc.extend(N[idx[done]].tolist()); idx = idx[~done]
            if idx.size == 0: break
        hit = s[N[idx] - p]
        best[idx[hit]] = p
        idx = idx[~hit]
    exc.extend(N[idx].tolist())               # uncovered after all starters < NMAX
    exc = sorted(exc)
    cov = best > 0
    imax = int(np.argmax(best))
    rows.append({
        "g": g,
        "odd_primes_dividing_g": "x".join(str(l) for l in range(3, g + 1, 2)
                                          if g % l == 0 and all(l % d for d in range(3, int(l**.5) + 1, 2))),
        "starters_below_NMAX": len(st),
        "density_vs_twin_measured": round(len(st) / twin_count, 4),
        "density_vs_twin_HL": round(hl_factor(g), 4),
        "total_exceptions": len(exc),
        "last_exception": exc[-1] if exc else "",
        "exceptions": " ".join(map(str, exc)) if len(exc) <= 25 else "(many)",
        "mean_min_starter": round(float(best[cov].mean()), 2),
        "max_min_starter": int(best[imax]),
        "N_at_max_min_starter": int(N[imax]),
        "seconds": round(time.time() - t0, 1),
    })
    r = rows[-1]
    print(f"g={g:3d} exc={r['total_exceptions']:3d} last={r['last_exception']!s:>6} "
          f"mean={r['mean_min_starter']:7.2f} max={r['max_min_starter']:6d} "
          f"dens={r['density_vs_twin_measured']:.3f}/{r['density_vs_twin_HL']:.3f} t={r['seconds']}s",
          flush=True)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"twin starters below {NMAX}: {twin_count}; wrote {OUT}")
