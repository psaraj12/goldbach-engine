"""Model C': Li2(M_g(X)) = C * (ln X)^2 / S(g), with Li2(y) = int_2^y dt / ln^2 t.

Only the constant C is fitted (on X = 1e4..1e6); the exponents are the heuristic
values 2 and 1. Also fits the free-exponent version B': ln Li2(M) = a + alpha lnln X - beta ln S.
Usage: python3 fit_li2.py gap_growth_1e8.csv
"""
import sys, csv
import numpy as np
from scipy.integrate import quad

R = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "gap_growth_1e8.csv")))
X = np.array([float(r["X"]) for r in R]); S = np.array([float(r["hl_factor"]) for r in R])
M = np.array([float(r["M"]) for r in R])

grid = np.unique(np.concatenate([np.arange(3, 20001), np.geomspace(20001, 2e6, 4000)]))
li2 = np.concatenate([[0], np.cumsum([quad(lambda t: 1 / np.log(t) ** 2, a, b)[0]
                                      for a, b in zip(grid[:-1], grid[1:])])])
Li2 = lambda y: np.interp(y, grid, li2)
Li2inv = lambda v: np.interp(v, li2, grid)

fit = np.isin(X, [1e4, 1e5, 1e6])
lnC = np.mean(np.log(Li2(M[fit])) - 2 * np.log(np.log(X[fit])) + np.log(S[fit]))
pred = Li2inv(np.exp(lnC) * np.log(X) ** 2 / S)
print(f"Model C': Li2(M) = C (ln X)^2 / S,  C = {np.exp(lnC):.4f}")
for L in sorted(set(X)):
    m = X == L; e = (pred[m] - M[m]) / M[m]
    tag = "fit" if L in (1e4, 1e5, 1e6) else ("held out" if L > 1e6 else "not fitted")
    print(f"  X=1e{int(round(np.log10(L)))} ({tag}): mean error {100*e.mean():+5.1f}%, "
          f"median |error| {100*np.median(abs(e)):4.1f}%")

D = np.column_stack([np.ones(fit.sum()), np.log(np.log(X[fit])), -np.log(S[fit])])
c = np.linalg.lstsq(D, np.log(Li2(M[fit])), rcond=None)[0]
print(f"Model B' (free exponents): alpha = {c[1]:.2f}, beta = {c[2]:.2f}")
