"""Fit growth laws for the largest minimal g-starter M_g(X).

Model A (proposed):  ln M = a + alpha*lnln X - beta*ln S(g)
Model B (heuristic): ln M - 2 lnln M = a + alpha*lnln X - beta*ln S(g)
   (from: P(no starter <= y works for N) ~ exp(-c S y / (ln^2 y ln N)),
    extreme over ~X/2 values of N  =>  S*M/ln^2 M ~ ln^2 X)
Fit on X in FIT_X, report errors on held-out limits.
"""
import sys, csv
import numpy as np

src = sys.argv[1]
rows = list(csv.DictReader(open(src)))
g = np.array([int(r["g"]) for r in rows]); X = np.array([float(r["X"]) for r in rows])
S = np.array([float(r["hl_factor"]) for r in rows]); M = np.array([float(r["M"]) for r in rows])
limits = sorted(set(X))
FIT_X = [1e4, 1e5, 1e6]

def design(m):
    return np.column_stack([np.ones(m.sum()), np.log(np.log(X[m])), -np.log(S[m])])

for name, ytrans, yinv in [
    ("A: ln M = a + alpha lnln X - beta ln S", lambda M: np.log(M), None),
    ("B: ln(M/ln^2 M) = a + alpha lnln X - beta ln S", lambda M: np.log(M) - 2*np.log(np.log(M)), "B")]:
    fit = np.isin(X, FIT_X)
    coef, *_ = np.linalg.lstsq(design(fit), ytrans(M[fit]), rcond=None)
    a, alpha, beta = coef
    yhat = design(np.ones_like(fit)) @ coef
    if yinv == "B":   # invert ln M - 2 lnln M = y by fixed-point iteration
        lm = yhat.copy()
        for _ in range(60): lm = yhat + 2*np.log(lm)
        pred = np.exp(lm)
    else:
        pred = np.exp(yhat)
    resid_fit = ytrans(M[fit]) - (design(fit) @ coef)
    r2 = 1 - resid_fit.var() / ytrans(M[fit]).var()
    print(f"\nModel {name}\n  fitted on X = 1e4..1e6: alpha = {alpha:.2f}, beta = {beta:.2f}, R^2 = {r2:.3f}")
    for L in limits:
        m = X == L
        err = (pred[m] - M[m]) / M[m]
        tag = "fit" if L in FIT_X else ("held out" if L > max(FIT_X) else "not fitted")
        print(f"  X = 1e{int(round(np.log10(L)))} ({tag:10s}): median |error| = {100*np.median(np.abs(err)):5.1f}%, "
              f"mean error = {100*err.mean():+6.1f}%, actual mean M = {M[m].mean():7.0f}, predicted mean M = {pred[m].mean():7.0f}")

# dispersion across gaps vs density factor, at the largest X
Lmax = max(limits); m = X == Lmax
print(f"\nSpread of M_g({Lmax:.0e}) across gaps, by density factor (log scale removes the size effect):")
print("  S group     | gaps | mean M | sd of ln M")
groups = [("2.0", lambda s: abs(s-2.0) < 1e-6), ("2.02-2.2", lambda s: 2.0 < s < 2.2),
          ("2.2-2.5", lambda s: 2.2 <= s < 2.5), ("2.6-2.7", lambda s: 2.6 <= s < 2.7), (">2.7", lambda s: s >= 2.7)]
for lab, f in groups:
    sel = m & np.array([f(v) for v in S])
    if sel.sum() > 1:
        print(f"  {lab:11s} | {sel.sum():4d} | {M[sel].mean():6.0f} | {np.log(M[sel]).std(ddof=1):.3f}")
