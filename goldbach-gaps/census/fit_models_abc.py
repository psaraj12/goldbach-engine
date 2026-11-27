"""Models A, B, C for M_g(X); fit on X = 1e4..1e6, test on 1e7, 1e8.
A: ln M          = a + alpha lnln X - beta ln S   (alpha, beta free)
B: ln M - 2lnlnM = a + alpha lnln X - beta ln S   (alpha, beta free)
C: ln M - 2lnlnM = a + 2 lnln X - 1 ln S          (only a fitted)
95% intervals for alpha, beta: bootstrap over gaps (2000 resamples), which
keeps each gap's nested M_g(X) values together."""
import csv, numpy as np
R = list(csv.DictReader(open("gap_growth_1e8.csv")))
G = np.array([int(r["g"]) for r in R]); X = np.array([float(r["X"]) for r in R])
S = np.array([float(r["hl_factor"]) for r in R]); M = np.array([float(r["M"]) for r in R])
FIT = [1e4, 1e5, 1e6]; TEST = [1e7, 1e8]
lx, ls = np.log(np.log(X)), np.log(S)
yA, yB = np.log(M), np.log(M) - 2*np.log(np.log(M))

def fit(y, m, model):
    if model == "C":
        a = np.mean(y[m] - 2*lx[m] + ls[m]); return np.array([a, 2.0, 1.0])
    D = np.column_stack([np.ones(m.sum()), lx[m], -ls[m]])
    return np.linalg.lstsq(D, y[m], rcond=None)[0]

def predict(coef, model):
    yhat = coef[0] + coef[1]*lx - coef[2]*ls
    if model == "A": return np.exp(yhat)
    lm = yhat.copy()
    for _ in range(80): lm = yhat + 2*np.log(lm)
    return np.exp(lm)

fitmask = np.isin(X, FIT)
rng = np.random.default_rng(1); gaps = np.unique(G)
for model, y in (("A", yA), ("B", yB), ("C", yB)):
    c = fit(y, fitmask, model); p = predict(c, model)
    line = f"Model {model}: alpha={c[1]:.2f} beta={c[2]:.2f}"
    if model != "C":
        boots = []
        for _ in range(2000):
            pick = rng.choice(gaps, len(gaps), replace=True)
            idx = np.concatenate([np.nonzero((G == g) & fitmask)[0] for g in pick])
            m = np.zeros(len(G), bool); m[idx] = True
            D = np.column_stack([np.ones(len(idx)), lx[idx], -ls[idx]])
            boots.append(np.linalg.lstsq(D, y[idx], rcond=None)[0])
        b = np.array(boots)
        line += (f"  95% CI alpha [{np.percentile(b[:,1],2.5):.2f}, {np.percentile(b[:,1],97.5):.2f}]"
                 f"  beta [{np.percentile(b[:,2],2.5):.2f}, {np.percentile(b[:,2],97.5):.2f}]")
    print(line)
    for L in TEST:
        m = X == L; e = (p[m] - M[m]) / M[m]
        print(f"   X=1e{int(np.log10(L))}: mean error {100*e.mean():+5.1f}%, median |error| {100*np.median(abs(e)):4.1f}%")
