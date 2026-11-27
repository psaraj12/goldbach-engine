#!/usr/bin/env python3
"""Pre-registered predictions for the minimal-starter window experiment.

Statistic: M = largest minimal g-starter partner over a window of W even N at
height N (window start; the window is tiny relative to N).

Four candidate laws, each "size(M) = C * count", calibrated ONCE on the
exhaustive census (gap_growth_1e8.csv, all X = 1e4..1e8, N = X, W = #even N <= X)
and then frozen. Nothing is refitted on window data.
  size:   K  = K_g(M), the exact number of g-starters <= M (sieve)
          L  = S(g) * Li2(M),  Li2(y) = int_2^y dt / ln^2 t
  count:  NW = ln N * ln W     (heuristic: q-primality ~ 1/ln N, extreme over W values)
          NN = (ln N)^2        (no dependence on window size; null for the W effect)
Census finding that motivates testing all four: with C fitted on 1e4..1e6, the
exact-count law K-NW drifts by +17% from 1e4 to 1e8 (the starters needed grow
faster than ln X ln W), while L-NW drifts by -5%; L partly cancels the drift
through its own over-count of small starters. The windows separate these.

Usage
  (from goldbach-gaps/)
  python3 prereg/prereg_windows.py predict --census census/gap_growth_1e8.csv --threads 12 --out prereg [--with-1e8]
      -> prereg/predictions.csv, prereg/run_grid.sh, prereg/PREREG.md (with SHA-256)
  python3 prereg/prereg_windows.py compare --pred prereg/predictions.csv --stats-dir runs
"""
import argparse, csv, glob, hashlib, math, os, sys, datetime
import numpy as np

GAPS    = [6, 30, 210, 2310, 30030, 510510]
HEIGHTS = [10**9, 10**12, 10**15, 10**18, 4 * 10**18]
TIERS   = [(10**5, 8), (10**6, 4), (10**7, 2)]     # (window size W in even N, replicates)
OPTIONAL_TIER = (10**8, 1)                         # long tier, --with-1e8
ANCHOR_LIMIT = 10**6
BAND = 0.80                                        # pre-registered prediction band coverage

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

def cumulative_integral(g, ymax=5e7, n=400000):
    """Grid t and cumulative int_2^t f, f = 1/(ln t ln(t+g)) (g=None: 1/ln^2 t)."""
    t = np.geomspace(2.0, ymax, n)
    f = 1.0 / (np.log(t) * np.log(t + g)) if g is not None else 1.0 / np.log(t) ** 2
    F = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(t))])
    return t, F

_cache = {}
def integ(g):
    if g not in _cache: _cache[g] = cumulative_integral(g)
    return _cache[g]
def I(g, y):     t, F = integ(g);    return np.interp(y, t, F)
def Iinv(g, v):  t, F = integ(g);    return np.interp(v, F, t)

MODELS = ["K-NW", "L-NW", "K-NN", "L-NN"]
PRIMARY = "K-NW"

_kcache = {}
def starter_counts(g, ymax=3 * 10**6):
    if g not in _kcache:
        top = ymax + g + 2
        sv = np.ones(top + 1, dtype=bool); sv[:2] = False
        for i in range(2, int(top ** 0.5) + 1):
            if sv[i]: sv[i * i::i] = False
        pr = np.nonzero(sv[:ymax])[0]; pr = pr[pr > 2]
        st = pr[sv[pr + g]]
        _kcache[g] = st
    return _kcache[g]
def K(g, y):      return float(np.searchsorted(starter_counts(g), y, side="right"))
def Kinv(g, v):
    st = starter_counts(g); i = int(math.ceil(v)) - 1
    return float(st[min(max(i, 0), len(st) - 1)])
def size(model, g, S, M):  return K(g, M) if model[0] == "K" else S * I(None, M)
def size_inv(model, g, S, v): return Kinv(g, v) if model[0] == "K" else Iinv(None, v / S)
def count(model, lnN, lnW): return lnN * lnW if model.endswith("NW") else lnN * lnN

def calibrate(census_csv):
    rows = [(int(r["g"]), float(r["X"]), float(r["M"]), float(r["hl_factor"]))
            for r in csv.DictReader(open(census_csv)) if float(r["X"]) >= 1e4]
    out = {}
    for m in MODELS:
        def lr(g, X, M, S, C):
            pred = size_inv(m, g, S, C * count(m, math.log(X), math.log((X - 6) / 2 + 1)))
            return math.log(M / pred)
        lnC = np.mean([math.log(size(m, g, S, M)) - math.log(count(m, math.log(X), math.log((X - 6) / 2 + 1)))
                       for g, X, M, S in rows])
        C = math.exp(lnC)
        res = [lr(g, X, M, S, C) for g, X, M, S in rows]
        byX = {k: math.exp(np.mean([lr(g, X, M, S, C) for g, X, M, S in rows if X == 10**k])) for k in range(4, 9)}
        q = np.quantile(res, [(1 - BAND) / 2, 1 - (1 - BAND) / 2])
        out[m] = dict(C=C, byX=byX, band=(float(q[0]), float(q[1])))
    return out

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

def predict(a):
    global GAPS, HEIGHTS, TIERS
    if a.smoke:   # pipeline test on windows OUTSIDE the registered grid
        GAPS, HEIGHTS, TIERS = [6, 2310], [10**10, 10**13], [(10**4, 2)]
    cal = calibrate(a.census)
    os.makedirs(a.out, exist_ok=True)
    tiers = TIERS + ([OPTIONAL_TIER] if a.with_1e8 else [])
    pred_csv = os.path.join(a.out, "predictions.csv")
    run_sh = os.path.join(a.out, "run_grid.sh")
    rows = []
    for g in GAPS:
        S = hl_factor(g)
        for H in HEIGHTS:
            start = H
            for W, reps in tiers:
                for rep in range(reps):
                    lo_N, hi_N = start, start + 2 * W - 2
                    lnN, lnW = 0.5 * (math.log(lo_N) + math.log(hi_N)), math.log(W)
                    bb = max(10, min(24, int(math.floor(math.log2(W / (4 * a.threads))))))
                    r = dict(gap=g, S=round(S, 6), height=H, W=W, rep=rep,
                             range_start=lo_N, range_end=hi_N, block_bits=bb)
                    for m in MODELS:
                        pm = size_inv(m, g, S, cal[m]["C"] * count(m, lnN, lnW))
                        bl, bh = cal[m]["band"]
                        r[f"pred_{m}"] = round(pm, 1)
                        r[f"lo_{m}"] = round(pm * math.exp(bl), 1)
                        r[f"hi_{m}"] = round(pm * math.exp(bh), 1)
                    rows.append(r)
                    start = hi_N + 2          # replicates: contiguous, disjoint
    with open(pred_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    with open(run_sh, "w") as f:
        f.write("#!/usr/bin/env bash\n# Generated by prereg_windows.py -- runs every pre-registered window.\n"
                "# Run from goldbach-gaps/:   ENGINE_DIR=engine THREADS=12 prereg/run_grid.sh\n"
                "# Results go to ./runs/. Windows already present in a minpstats file are skipped,\n"
                "# so the script can be re-run after an interruption.\n"
                "# Build once per gap, from goldbach-gaps/ (plain build is fine: the minimal-partner scan dominates):\n")
        for g in GAPS:
            f.write(f"#   g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX2 -DGAP={g} "
                    f"engine/goldbach_v6_4_gap.cpp -o engine/gb_gap{g}\n")
        f.write(f"set -euo pipefail\nTHREADS=${{THREADS:-{a.threads}}}\n"
                "ENGINE_DIR=$(cd \"${ENGINE_DIR:-engine}\" && pwd)\nmkdir -p runs && cd runs\n")
        for r in rows:
            f.write(f"[ -n \"$(awk -F, '$2=={r['range_start']} && $3=={r['range_end']}' "
                    f"goldbach_gap{r['gap']}_minpstats.csv 2>/dev/null)\" ] || "
                    f"\"$ENGINE_DIR/gb_gap{r['gap']}\" {r['range_start']} {r['range_end']} "
                    f"{r['block_bits']} $THREADS {ANCHOR_LIMIT} 0 1024 "
                    f"--minp-per-window={r['W']} --minp-limit=0 "
                    f"> log_gap{r['gap']}_{r['range_start']}.txt 2>&1\n")
    os.chmod(run_sh, 0o755)

    digest = sha256(pred_csv)
    n_even = sum(r["W"] for r in rows)
    drift = "\n".join(f"| {m} | {cal[m]['C']:.6g} | " + " | ".join(f"{cal[m]['byX'][k]:.3f}" for k in range(4, 9)) + " |"
                      for m in MODELS)
    with open(os.path.join(a.out, "PREREG.md"), "w") as f:
        f.write(f"""# Pre-registration: minimal-starter window experiment

Written: {datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
predictions.csv SHA-256: `{digest}`
Census used for calibration: `{os.path.basename(a.census)}` SHA-256 `{sha256(a.census)}`

Commit this file and predictions.csv BEFORE running run_grid.sh.

## Design
- Gaps: {GAPS}
- Heights (window start): {', '.join(f'{h:.0e}' for h in HEIGHTS)}
- Window sizes W (even N) x replicates: {', '.join(f'{w:.0e} x {k}' for w, k in tiers)}
- Windows: {len(rows)}; even N examined: {n_even:,}
- Engine: engine/goldbach_v6_4_gap.cpp built with -DGAP=g; p_anchor_limit = {ANCHOR_LIMIT}; statistic = max minimal g-starter partner over the window (one minpstats row per window).

## Frozen laws: size(M) = C x count, C calibrated on the census (X = 1e4..1e8), never refitted
- size K = exact number of g-starters <= M; size L = S(g) Li2(M)
- count NW = ln N ln W (heuristic); count NN = (ln N)^2 (null: no window-size effect)
- Primary law: {PRIMARY} (the heuristic with no approximation in the starter count)

Census residuals (geometric mean of actual/predicted M, by X), constants as frozen:

| law | C | 1e4 | 1e5 | 1e6 | 1e7 | 1e8 |
|---|---|---|---|---|---|---|
{drift}

{int(BAND*100)}% bands come from each law's census residuals.

## Hypotheses and success criteria (decided now, before any window is run)
1. **Window size matters (H1):** each NW law has a smaller mean |log error| than its NN partner (K-NW vs K-NN, L-NW vs L-NN).
2. **Drift continues (H2):** for the primary law K-NW, the geometric-mean actual/predicted rises with height: the value at 4e18 exceeds the value at 1e9.
3. **Accuracy:** report each law's geometric-mean actual/predicted, mean |log error| and band coverage; a law counts as adequate if its overall ratio lies in 0.85-1.15 and 65-95% of windows fall in its band.
4. **Gap uniformity:** for the best law, no gap has a geometric-mean ratio outside 0.80-1.25.

Each hypothesis is reported as SUPPORTED or NOT SUPPORTED as stated, whatever the outcome.
Any later model is reported separately as exploratory, not as a test.
""")
    print(f"wrote {pred_csv} ({len(rows)} windows, {n_even:,} even N), {run_sh}, {os.path.join(a.out, 'PREREG.md')}")
    for m in MODELS:
        print(f"  {m}: C = {cal[m]['C']:.6g}  census actual/pred by X: " +
              " ".join(f"{cal[m]['byX'][k]:.3f}" for k in range(4, 9)))
    print(f"predictions.csv SHA-256 {digest}")

def compare(a):
    preds = {(int(r["gap"]), int(r["range_start"]), int(r["range_end"])): r
             for r in csv.DictReader(open(a.pred))}
    got = {}
    for path in glob.glob(os.path.join(a.stats_dir, "goldbach_gap*_minpstats.csv")):
        for r in csv.DictReader(open(path)):
            g = int(r["prog_version"].split("gap")[1].split("_")[0])
            if r["resumed"] == "1" or int(r["n_without_partner"]) != 0: continue
            got[(g, int(r["range_start"]), int(r["range_end"]))] = float(r["max_min_p"])
    res = []
    for k, p in preds.items():
        if k not in got: continue
        M = got[k]
        e = dict(g=k[0], H=int(p["height"]), W=int(p["W"]))
        for m in MODELS:
            e[m] = math.log(M / float(p[f"pred_{m}"]))
            e["in_" + m] = float(p[f"lo_{m}"]) <= M <= float(p[f"hi_{m}"])
        res.append(e)
    if not res: sys.exit("no completed windows matched the predictions")
    gm = lambda xs: math.exp(np.mean(xs)) if xs else float("nan")
    print(f"windows scored: {len(res)} of {len(preds)}\n")
    print(f"{'law':6} {'actual/pred':>11} {'mean|log err|':>14} {'in band':>8}  adequate?")
    stat = {}
    for m in MODELS:
        r_ = gm([e[m] for e in res]); ae = np.mean([abs(e[m]) for e in res]); cov = np.mean([e["in_" + m] for e in res])
        stat[m] = ae
        ok = 0.85 <= r_ <= 1.15 and 0.65 <= cov <= 0.95
        print(f"{m:6} {r_:11.3f} {ae:14.3f} {100*cov:7.1f}%  {'yes' if ok else 'no'}")
    h1 = stat["K-NW"] < stat["K-NN"] and stat["L-NW"] < stat["L-NN"]
    print(f"\nH1 window size matters: {'SUPPORTED' if h1 else 'NOT SUPPORTED'}")
    Hs = sorted({e["H"] for e in res})
    lo = [e[PRIMARY] for e in res if e["H"] == Hs[0]]; hi = [e[PRIMARY] for e in res if e["H"] == Hs[-1]]
    if lo and hi:
        print(f"H2 drift continues ({PRIMARY}): ratio at lowest height {gm(lo):.3f}, at highest {gm(hi):.3f} -> "
              f"{'SUPPORTED' if gm(hi) > gm(lo) else 'NOT SUPPORTED'}")
    best = min(MODELS, key=lambda m: stat[m])
    print(f"\nbest law by mean |log error|: {best}")
    for key in ("g", "H", "W"):
        for v in sorted({e[key] for e in res}):
            print(f"   {key}={v:<22} actual/pred ({best}) = {gm([e[best] for e in res if e[key] == v]):.3f}")
    worst = [0.80 <= gm([e[best] for e in res if e["g"] == g]) <= 1.25 for g in sorted({e["g"] for e in res})]
    print(f"gap uniformity ({best}): {'PASS' if all(worst) else 'FAIL'}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("predict"); p.add_argument("--census", required=True)
    p.add_argument("--threads", type=int, default=12); p.add_argument("--out", default="prereg")
    p.add_argument("--with-1e8", action="store_true", help="add the long tier W = 1e8 (one window per gap and height)")
    p.add_argument("--smoke", action="store_true", help="tiny test grid at heights outside the registered grid (pipeline check only)")
    c = sub.add_parser("compare"); c.add_argument("--pred", required=True); c.add_argument("--stats-dir", required=True)
    a = ap.parse_args(); predict(a) if a.cmd == "predict" else compare(a)
