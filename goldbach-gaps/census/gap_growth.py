"""Growth of the largest minimal g-starter with the verification limit X.

For every gap g in GAPS, compute the least g-starter partner of every even N in
[6, XMAX] (exhaustive, no anchor limit), then report, for each X in LIMITS:
  M_g(X)       = max over representable even N <= X of the least starter
  N_at_M       = the N attaining it (smallest such N)
  mean_min     = mean least starter over representable even N <= X
  exceptions   = number of even N <= X with no g-starter partition
Output: long-format CSV, one row per (g, X).
"""
import sys, csv, time
import numpy as np

XMAX   = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
OUT    = sys.argv[2] if len(sys.argv) > 2 else "gap_growth.csv"
GAPS   = range(6, 601, 6)
LIMITS = [10**k for k in range(3, 9) if 10**k <= XMAX]

top = XMAX + max(GAPS) + 2
s = np.ones(top + 1, dtype=bool); s[:2] = False
for i in range(2, int(top ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
cand = np.nonzero(s[3:XMAX])[0].astype(np.int64) + 3          # odd primes < XMAX

nN = (XMAX - 6) // 2 + 1                                      # N = 6 + 2*i, i < nN
cut = {X: (X - 6) // 2 + 1 for X in LIMITS}                   # number of N <= X

def hl_factor(g):
    f, m, l = 1.0, g, 3
    while m % 2 == 0: m //= 2
    while l * l <= m:
        if m % l == 0:
            f *= (l - 1) / (l - 2)
            while m % l == 0: m //= l
        l += 2
    if m > 1: f *= (m - 1) / (m - 2)
    return f

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["g", "hl_factor", "X", "M", "N_at_M", "mean_min", "exceptions"])
    for g in GAPS:
        t0 = time.time()
        st = cand[s[cand + g]]
        best = np.zeros(nN, dtype=np.int32)
        idx = np.arange(nN, dtype=np.int64)
        for p in st:
            if idx.size == 0: break
            Nv = 6 + 2 * idx
            keep = Nv > p                       # p must be < N
            if not keep.all():
                idx = idx[keep]; Nv = Nv[keep]
                if idx.size == 0: break
            hit = s[Nv - p]
            best[idx[hit]] = p
            idx = idx[~hit]
        h = hl_factor(g)
        for X in LIMITS:
            b = best[:cut[X]]
            i = int(np.argmax(b))
            cov = b > 0
            w.writerow([g, round(h, 6), X, int(b[i]), 6 + 2 * i,
                        round(float(b[cov].mean()), 3), int((~cov).sum())])
        fh.flush()
        print(f"g={g:3d} M(1e8)={int(best.max()):5d} t={time.time()-t0:.1f}s", flush=True)
print("done", OUT)
