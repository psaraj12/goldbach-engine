# Cross-check against T. Pandey, arXiv:2608.02381v2, Section 6 (record table).
# Usage: python3 crosscheck_pandey_records.py   (sympy; runs in about a second)
from sympy import isprime
# Pandey's 51 record pairs (n, t_min), arXiv:2608.02381v2, Section 6
rec = [(6,3),(12,5),(30,7),(98,19),(220,29),(308,31),(556,107),(796,137),(2776,197),(6908,199),(8138,271),
(14318,1021),(66818,1231),(137528,1291),(344018,1321),(479186,1609),(528968,1621),(774968,1669),(1142224,1871),
(2214848,2269),(6137438,2689),(8753956,2999),(15704356,3299),(24968056,3359),(31480496,3469),(36612068,3919),
(59892398,4051),(65712758,4639),(241609594,5417),(453267964,5477),(639267856,5849),(1200364448,6271),
(1311113198,6451),(3104283256,7349),(3760188526,7457),(6082985476,7559),(8126254564,10271),(68607056186,11059),
(77739196118,14449),(571714791706,14549),(1821067383866,14629),(1855623735874,15647),(2435159540516,16189),
(3824295471338,16651),(4968233746816,17207),(5276047389878,17419),(6053540148034,18311),(6527240154856,20747),
(15149851353116,21493),(29330655440536,21599),(31819171719758,23029)]
twinmem = lambda t: isprime(t) and (isprime(t-2) or isprime(t+2))
starters = [p for p in range(3,40000) if isprime(p) and (isprime(p+2) or isprime(p+4))]
def least_tc(n): return next(p for p in starters if p < n and isprime(n-p))
def least_twinmem(n): return next(t for t in range(3,40000) if twinmem(t) and isprime(n-t))
agree=disagree=0; rows=[]
for n,t in rec:
    assert least_twinmem(n)==t, ("Pandey record not reproduced", n)
    s = least_tc(n)
    if n % 3 == 1:
        ok = (s == t); agree += ok; disagree += (not ok)
    rows.append((n, n%3, t, s))
print(f"all 51 of Pandey's records reproduced (least twin member, minimality checked)")
print(f"class N=1 (mod 3): {agree} records, least TC starter == his t_min in all? {disagree==0}")
print(f"{'n':>18} {'n%3':>3} {'his t_min':>9} {'our least TC starter':>20}")
for n,c,t,s in rows:
    if n > 10**7: print(f"{n:>18,} {c:>3} {t:>9,} {s:>20,}{'   (same: class 1)' if c==1 else ''}")
