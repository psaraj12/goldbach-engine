#!/usr/bin/env python3
"""Round 2 pre-registration: minimal-starter windows with the worst-N residue correction.

Statistic per window: M = largest minimal g-starter partner over W even N at height N;
also the N attaining it (the "hardest N").

Frozen laws (constants fitted ONCE on the 420 Round-1 windows, then never refitted):
  R2-F  (primary): K_g(M) = C * ln N * ln W / F(g)
  R2-FH          : K_g(M) = C * (ln N)^a * ln W / F(g)      (a fitted on Round 1)
  K-NW  (baseline): K_g(M) = C * ln N * ln W                  (Round-1 law, recalibrated on Round 1)
  K_g(M) = exact number of g-starters <= M (sieve)
  F(g)   = prod over primes 5 <= l <= 10^5, l not dividing g, of (1 - 1/(l-2)^2)
Mechanism statistics (no fitting):
  unfavourable(N, g, l): N != 0 (mod l), and also N != -g (mod l) when l does not divide g
  H_g(N) = prod over odd primes l <= 97 of (1 - b/l) / ((1 - a/l)(1 - 1/l)),
           a = |{0, -g} mod l|, b = |{0, -g, N} mod l|
           (three-form singular series of p, p+g, N-p over the pair series of p, p+g)

Usage (from goldbach-gaps/)
  python3 prereg_r2/prereg_r2.py predict --r1-pred prereg/predictions.csv --r1-runs runs --threads 12 --out prereg_r2
  python3 prereg_r2/prereg_r2.py compare --pred prereg_r2/predictions_r2.csv --stats-dir runs_r2
"""
import argparse, csv, glob, hashlib, math, os, sys, datetime
import numpy as np

GAPS    = [6, 12, 42, 66, 330, 390, 2310, 4620]
HEIGHTS = [3 * 10**9, 3 * 10**12, 3 * 10**15, 2 * 10**18, 8 * 10**18]
TIERS   = [(10**5, 8), (10**6, 4), (10**7, 2)]
ANCHOR_LIMIT = 10**6
BAND = 0.80
LAWS = ["R2-F", "R2-FH", "K-NW"]
PRIMARY = "R2-F"
H_PRIMES_MAX = 97          # truncation of H_g(N)
H_SAMPLE = 1000            # random even N per window for H5
H4_PRIMES = [3, 5, 7, 11, 13]
ALPHA = 0.001

# ---------------- arithmetic helpers ----------------
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR_F = primes_upto(10**5)
PR_H = [int(l) for l in primes_upto(H_PRIMES_MAX) if l >= 3]

def F(g):
    return float(np.prod([1 - 1 / (l - 2) ** 2 for l in PR_F if l >= 5 and g % l]))

_kc = {}
def starters(g, ymax=3 * 10**6):
    if g not in _kc:
        top = ymax + g + 2
        sv = np.ones(top + 1, dtype=bool); sv[:2] = False
        for i in range(2, int(top ** 0.5) + 1):
            if sv[i]: sv[i * i::i] = False
        pr = np.nonzero(sv[:ymax])[0]; pr = pr[pr > 2]
        _kc[g] = pr[sv[pr + g]]
    return _kc[g]
def K(g, y):    return float(np.searchsorted(starters(g), y, side="right"))
def Kinv(g, v):
    st = starters(g); i = int(math.ceil(v)) - 1
    return float(st[min(max(i, 0), len(st) - 1)])

def unfavourable(N, g, l):
    if g % l == 0: return N % l != 0
    return N % l != 0 and (N + g) % l != 0

def H(N, g):
    h = 1.0
    for l in PR_H:
        a = len({0, (-g) % l}); b = len({0, (-g) % l, N % l})
        h *= (1 - b / l) / ((1 - a / l) * (1 - 1 / l))
    return h

def window_rng(g, lo, hi):
    seed = int(hashlib.sha256(f"r2:{g}:{lo}:{hi}".encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

# ---------------- reading results ----------------
def read_stats(stats_dir):
    out = {}
    for path in sorted(glob.glob(os.path.join(stats_dir, "goldbach_gap*_minpstats.csv"))):
        for r in csv.DictReader(open(path)):
            g = int(r["prog_version"].split("gap")[1].split("_")[0])
            if r["resumed"] == "1" or int(r["n_without_partner"]) != 0: continue
            out[(g, int(r["range_start"]), int(r["range_end"]))] = (float(r["max_min_p"]), int(r["N_at_max"]))
    return out

def law_value(law, g, lnN, lnW, P):
    if law == "K-NW":  return P["K-NW"]["C"] * lnN * lnW
    if law == "R2-F":  return P["R2-F"]["C"] * lnN * lnW / F(g)
    return P["R2-FH"]["C"] * lnN ** P["R2-FH"]["a"] * lnW / F(g)

# ---------------- calibration on Round 1 ----------------
def calibrate(r1_pred, r1_runs):
    meta = {(int(r["gap"]), int(r["range_start"]), int(r["range_end"])): r for r in csv.DictReader(open(r1_pred))}
    got = read_stats(r1_runs)
    data = []
    for k, (M, Nmax) in got.items():
        if k not in meta: continue
        r = meta[k]
        lnN = 0.5 * (math.log(k[1]) + math.log(k[2])); lnW = math.log(int(r["W"]))
        data.append((k[0], lnN, lnW, M))
    if len(data) != 420: sys.exit(f"expected the 420 Round-1 windows, found {len(data)}")
    P = {}
    P["K-NW"] = dict(C=math.exp(np.mean([math.log(K(g, M)) - math.log(lnN * lnW) for g, lnN, lnW, M in data])))
    P["R2-F"] = dict(C=math.exp(np.mean([math.log(K(g, M) * F(g)) - math.log(lnN * lnW) for g, lnN, lnW, M in data])))
    y = np.array([math.log(K(g, M) * F(g) / lnW) for g, lnN, lnW, M in data])
    X = np.column_stack([np.ones(len(data)), [math.log(lnN) for g, lnN, lnW, M in data]])
    c, a = np.linalg.lstsq(X, y, rcond=None)[0]
    P["R2-FH"] = dict(C=math.exp(c), a=float(a))
    for law in LAWS:
        res = [math.log(M / Kinv(g, law_value(law, g, lnN, lnW, P))) for g, lnN, lnW, M in data]
        q = np.quantile(res, [(1 - BAND) / 2, 1 - (1 - BAND) / 2])
        P[law]["band"] = (float(q[0]), float(q[1]))
        P[law]["r1_ratio"] = math.exp(np.mean(res)); P[law]["r1_abs"] = float(np.mean(np.abs(res)))
    return P

# ---------------- predict ----------------
def predict(a):
    global GAPS, HEIGHTS, TIERS
    if a.smoke:
        GAPS, HEIGHTS, TIERS = [12, 4620], [5 * 10**10, 5 * 10**13], [(10**4, 2)]
    P = calibrate(a.r1_pred, a.r1_runs)
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for g in GAPS:
        for Hh in HEIGHTS:
            start = Hh
            for W, reps in TIERS:
                for rep in range(reps):
                    lo, hi = start, start + 2 * W - 2
                    lnN, lnW = 0.5 * (math.log(lo) + math.log(hi)), math.log(W)
                    bb = max(10, min(24, int(math.floor(math.log2(W / (4 * a.threads))))))
                    r = dict(gap=g, F=round(F(g), 6), height=Hh, W=W, rep=rep,
                             range_start=lo, range_end=hi, block_bits=bb)
                    for law in LAWS:
                        pm = Kinv(g, law_value(law, g, lnN, lnW, P)); bl, bh = P[law]["band"]
                        r[f"pred_{law}"] = round(pm, 1)
                        r[f"lo_{law}"] = round(pm * math.exp(bl), 1); r[f"hi_{law}"] = round(pm * math.exp(bh), 1)
                    r["p_unfav_" + "_".join(map(str, H4_PRIMES))] = round(
                        float(np.prod([((l - 1) / l) if g % l == 0 else ((l - 2) / l) for l in H4_PRIMES])), 6)
                    rows.append(r)
                    start = hi + 2
    pred_csv = os.path.join(a.out, "predictions_r2.csv")
    with open(pred_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    run_sh = os.path.join(a.out, "run_grid_r2.sh")
    with open(run_sh, "w") as f:
        f.write("#!/usr/bin/env bash\n# Generated by prereg_r2.py -- runs every Round-2 window.\n"
                "# Run from goldbach-gaps/:   ENGINE_DIR=engine THREADS=12 prereg_r2/run_grid_r2.sh\n"
                "# Results go to ./runs_r2/ (kept apart from Round 1). Finished windows are skipped.\n"
                "# Build once per gap, from goldbach-gaps/:\n")
        for g in GAPS:
            f.write(f"#   g++ -O3 -march=native -fopenmp -std=c++17 -DSIMD_AVX2 -DGAP={g} "
                    f"engine/goldbach_v6_4_gap.cpp -o engine/gb_gap{g}\n")
        f.write(f"set -euo pipefail\nTHREADS=${{THREADS:-{a.threads}}}\n"
                "ENGINE_DIR=$(cd \"${ENGINE_DIR:-engine}\" && pwd)\nmkdir -p runs_r2 && cd runs_r2\n")
        for r in rows:
            f.write(f"[ -n \"$(awk -F, '$2=={r['range_start']} && $3=={r['range_end']}' "
                    f"goldbach_gap{r['gap']}_minpstats.csv 2>/dev/null)\" ] || "
                    f"\"$ENGINE_DIR/gb_gap{r['gap']}\" {r['range_start']} {r['range_end']} "
                    f"{r['block_bits']} $THREADS {ANCHOR_LIMIT} 0 1024 "
                    f"--minp-per-window={r['W']} --minp-limit=0 "
                    f"> log_gap{r['gap']}_{r['range_start']}.txt 2>&1\n")
    os.chmod(run_sh, 0o755)

    r1_files = sorted(glob.glob(os.path.join(a.r1_runs, "goldbach_gap*_minpstats.csv")))
    r1_hash = hashlib.sha256(b"".join(open(p, "rb").read() for p in r1_files)).hexdigest()
    digest = sha256(pred_csv); nN = sum(r["W"] for r in rows)
    lawtab = "\n".join(f"| {law} | {P[law]['C']:.6g} | {P[law].get('a', 1.0):.4f} | {P[law]['r1_ratio']:.3f} | {P[law]['r1_abs']:.3f} |" for law in LAWS)
    with open(os.path.join(a.out, "PREREG_R2.md"), "w") as f:
        f.write(f"""# Pre-registration, Round 2: worst-N residue correction for extreme starter witnesses

Written: {datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
predictions_r2.csv SHA-256: `{digest}`
Calibration input: the 420 Round-1 windows, minpstats files in `{a.r1_runs}/` (combined SHA-256 `{r1_hash}`),
with window metadata from `{a.r1_pred}` (SHA-256 `{sha256(a.r1_pred)}`).

Commit this file and predictions_r2.csv, and push, BEFORE running run_grid_r2.sh.

## Design
- Gaps: {GAPS} (6 and 2310 carried over from Round 1; 12 and 4620 share their odd prime divisors)
- Heights (window start): {', '.join(f'{h:.0e}' for h in HEIGHTS)} (none used in Round 1)
- Window sizes W x replicates: {', '.join(f'{w:.0e} x {k}' for w, k in TIERS)}
- Windows: {len(rows)}; even N examined: {nN:,}
- Engine: engine/goldbach_v6_4_gap.cpp with -DGAP=g, p_anchor_limit = {ANCHOR_LIMIT}; results in runs_r2/.

## Frozen laws (fitted on Round 1, never refitted)
F(g) = product over primes 5 <= l <= 10^5 with l not dividing g of (1 - 1/(l-2)^2).

| law | C | ln N exponent | Round-1 actual/pred | Round-1 mean abs log error |
|---|---|---|---|---|
{lawtab}

Primary law: {PRIMARY}. {int(BAND*100)}% bands come from each law's Round-1 residuals.

## Mechanism statistics (frozen)
- Unfavourable class: N != 0 (mod l), and N != -g (mod l) when l does not divide g. Primes tested: {H4_PRIMES}.
- H_g(N) = product over odd primes l <= {H_PRIMES_MAX} of (1 - b/l)/((1 - a/l)(1 - 1/l)),
  a = |{{0, -g}} mod l|, b = |{{0, -g, N}} mod l|.
- For each window: the hardest N (minpstats N_at_max) is ranked against {H_SAMPLE} random even N from the same
  window, drawn with numpy default_rng seeded by the first 16 hex digits of SHA-256("r2:g:start:end").
  Mid-rank = (#sample with H < H_hard + 0.5 * #ties) / {H_SAMPLE}.

## Hypotheses and criteria (decided now; each reported as SUPPORTED / NOT SUPPORTED)
- **H1 (gap correction helps):** R2-F has a smaller mean |log error| than K-NW over all Round-2 windows.
- **H2 (height law holds):** for R2-FH, the geometric-mean actual/predicted at 8e18 divided by that at 3e9 lies in 0.90-1.10.
- **H3 (same divisors, same behaviour):** under R2-F, the geometric-mean actual/predicted of g=12 over g=6, and of g=4620 over g=2310, each lie in 0.90-1.10.
- **H4 (hardest N are unfavourable):** for each l in {H4_PRIMES}, the number of windows whose hardest N is unfavourable at l
  exceeds its expectation under random even N, one-sided normal approximation to the Poisson-binomial, p < {ALPHA}.
  Supported if it holds for every l.
- **H5 (hardest N have small H):** the mean mid-rank of H_g(hardest N) is below 0.40, and more than half the windows have
  mid-rank below 0.5 (one-sided sign test, p < {ALPHA}).
- **Accuracy (reported for each law):** adequate if overall actual/predicted lies in 0.85-1.15 and 65-95% of windows fall in its band.

No further laws are added after registration; anything else is reported as exploratory.
""")
    print(f"wrote {pred_csv} ({len(rows)} windows, {nN:,} even N), {run_sh}, PREREG_R2.md")
    for law in LAWS:
        print(f"  {law}: C = {P[law]['C']:.6g}, a = {P[law].get('a', 1.0):.4f}, Round-1 ratio {P[law]['r1_ratio']:.3f}, mean |log err| {P[law]['r1_abs']:.3f}")
    print(f"predictions_r2.csv SHA-256 {digest}")

# ---------------- compare ----------------
def compare(a):
    from math import erf, sqrt
    preds = {(int(r["gap"]), int(r["range_start"]), int(r["range_end"])): r for r in csv.DictReader(open(a.pred))}
    got = read_stats(a.stats_dir)
    res = []
    for k, p in preds.items():
        if k not in got: continue
        M, Nh = got[k]; g, lo, hi = k
        e = dict(g=g, H=int(p["height"]), W=int(p["W"]), lo=lo, hi=hi, Nh=Nh)
        for law in LAWS:
            e[law] = math.log(M / float(p[f"pred_{law}"]))
            e["in_" + law] = float(p[f"lo_{law}"]) <= M <= float(p[f"hi_{law}"])
        res.append(e)
    if not res: sys.exit("no completed windows matched the predictions")
    gm = lambda xs: math.exp(np.mean(xs)) if xs else float("nan")
    pnorm = lambda z: 0.5 * (1 - erf(z / sqrt(2)))     # upper tail
    print(f"windows scored: {len(res)} of {len(preds)}\n")
    print(f"{'law':6} {'actual/pred':>11} {'mean|log err|':>14} {'in band':>8}  adequate?")
    ae = {}
    for law in LAWS:
        r_ = gm([e[law] for e in res]); ae[law] = np.mean([abs(e[law]) for e in res]); cov = np.mean([e["in_" + law] for e in res])
        print(f"{law:6} {r_:11.3f} {ae[law]:14.3f} {100*cov:7.1f}%  {'yes' if 0.85 <= r_ <= 1.15 and 0.65 <= cov <= 0.95 else 'no'}")
    print(f"\nH1 gap correction helps: R2-F {ae['R2-F']:.3f} vs K-NW {ae['K-NW']:.3f} -> {'SUPPORTED' if ae['R2-F'] < ae['K-NW'] else 'NOT SUPPORTED'}")
    Hs = sorted({e["H"] for e in res})
    h2 = gm([e["R2-FH"] for e in res if e["H"] == Hs[-1]]) / gm([e["R2-FH"] for e in res if e["H"] == Hs[0]])
    print(f"H2 height law (R2-FH): top/bottom height ratio {h2:.3f} -> {'SUPPORTED' if 0.90 <= h2 <= 1.10 else 'NOT SUPPORTED'}")
    pr = {g: gm([e[PRIMARY] for e in res if e["g"] == g]) for g in sorted({e["g"] for e in res})}
    pairs = [(12, 6), (4620, 2310)]
    pv = [(x, y, pr[x] / pr[y]) for x, y in pairs if x in pr and y in pr]
    for x, y, v in pv: print(f"   g={x} over g={y}: {v:.3f}")
    print(f"H3 same divisors: {'SUPPORTED' if pv and all(0.90 <= v <= 1.10 for _, _, v in pv) else 'NOT SUPPORTED'}")
    h4 = True
    for l in H4_PRIMES:
        obs = sum(unfavourable(e["Nh"], e["g"], l) for e in res)
        ps = [((l - 1) / l) if e["g"] % l == 0 else ((l - 2) / l) for e in res]
        mu, var = sum(ps), sum(q * (1 - q) for q in ps)
        z = (obs - mu) / math.sqrt(var); pval = pnorm(z); h4 &= pval < ALPHA
        print(f"   l={l:2d}: hardest N unfavourable in {obs}/{len(res)} windows, expected {mu:.1f}, z = {z:+.1f}, p = {pval:.2g}")
    print(f"H4 hardest N unfavourable: {'SUPPORTED' if h4 else 'NOT SUPPORTED'}")
    ranks = []
    for e in res:
        rng = window_rng(e["g"], e["lo"], e["hi"]); W = (e["hi"] - e["lo"]) // 2 + 1
        sample = [e["lo"] + 2 * int(i) for i in rng.integers(0, W, size=H_SAMPLE)]
        hh = H(e["Nh"], e["g"]); hs = np.array([H(n, e["g"]) for n in sample])
        ranks.append((np.sum(hs < hh - 1e-12) + 0.5 * np.sum(np.abs(hs - hh) <= 1e-12)) / H_SAMPLE)
    ranks = np.array(ranks); below = int(np.sum(ranks < 0.5)); n = len(ranks)
    pz = pnorm((below - n / 2) / math.sqrt(n / 4))
    print(f"   mean mid-rank of H_g(hardest N) = {ranks.mean():.3f}; {below}/{n} windows below 0.5 (sign test p = {pz:.2g})")
    print(f"H5 hardest N have small H: {'SUPPORTED' if ranks.mean() < 0.40 and pz < ALPHA else 'NOT SUPPORTED'}")
    for key in ("g", "H", "W"):
        for v in sorted({e[key] for e in res}):
            print(f"   {key}={v:<22} actual/pred ({PRIMARY}) = {gm([e[PRIMARY] for e in res if e[key] == v]):.3f}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("predict")
    p.add_argument("--r1-pred", required=True); p.add_argument("--r1-runs", required=True)
    p.add_argument("--threads", type=int, default=12); p.add_argument("--out", default="prereg_r2")
    p.add_argument("--smoke", action="store_true", help="tiny grid OUTSIDE the registered one (pipeline check)")
    c = sub.add_parser("compare"); c.add_argument("--pred", required=True); c.add_argument("--stats-dir", required=True)
    a = ap.parse_args(); predict(a) if a.cmd == "predict" else compare(a)
