"""Actual g-starter counts below y vs the Hardy-Littlewood predictions
   S(g)*C2*2 ... written as  K_g(y) ~ 2*C2*S'(g) * integral,  where S'(g) = prod_{l|g, l odd} (l-1)/(l-2)
   and the integral is either Li2(y) = int dt/ln^2 t  or  I_g(y) = int dt/(ln t ln(t+g)).
   Also reports the effective starter density relative to g = 6 at each y."""
import numpy as np
from scipy.integrate import quad
C2 = 0.6601618158468696
Y = [10**3, 10**4, 10**5, 10**6, 10**7]
GAPS = [6, 30, 210, 2310, 30030, 510510]
top = max(Y) + max(GAPS) + 2
s = np.ones(top + 1, bool); s[:2] = False
for i in range(2, int(top**.5) + 1):
    if s[i]: s[i*i::i] = False
pr = np.nonzero(s)[0]; pr = pr[pr > 2]
def Sp(g):
    f, m, l = 1.0, g, 3
    while m % 2 == 0: m //= 2
    while l*l <= m:
        if m % l == 0:
            f *= (l-1)/(l-2)
            while m % l == 0: m //= l
        l += 2
    if m > 1: f *= (m-1)/(m-2)
    return f
def integ(f, y):  # integral from 2 to y, split for accuracy
    pts = np.unique(np.concatenate([[2], np.geomspace(2, y, 60)]))
    return sum(quad(f, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
print(f"{'g':>7} {'S(g)':>6} | " + " | ".join(f"y=1e{int(np.log10(y))}: actual / Li2 / I_g" for y in Y))
count6 = {}
for g in GAPS:
    cells = []
    for y in Y:
        c = int(np.count_nonzero(s[pr[pr < y] + g]))
        if g == 6: count6[y] = c
        li2 = 2*C2*Sp(g)*integ(lambda t: 1/np.log(t)**2, y)
        ig  = 2*C2*Sp(g)*integ(lambda t: 1/(np.log(t)*np.log(t+g)), y)
        cells.append(f"{c:7d} / {c/li2:5.2f} / {c/ig:5.2f}")
    print(f"{g:7d} {Sp(g)*1:6.3f} | " + " | ".join(cells))
print("\n(ratios = actual / predicted; 1.00 = perfect)")
print("\nEffective starter density relative to g = 6 (actual counts):")
for g in GAPS:
    print(f"  g={g:7d}: " + "  ".join(f"y=1e{int(np.log10(y))}: {np.count_nonzero(s[pr[pr<y]+g])/count6[y]:5.2f}" for y in Y))
