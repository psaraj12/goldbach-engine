"""Independent audit of a Conjecture-B (sexy-starter) campaign.
Usage: python3 audit_sx.py [results_dir]   (default: results/conjecture_B)
Rechecks every sampled witness and every logged minimal witness with sympy,
spot-checks minimality, and reports the per-class (N mod 3) statistics."""
import re, statistics, sys, os
from sympy import isprime
D = sys.argv[1] if len(sys.argv) > 1 else "results/conjecture_B"
# ---- checkpoint ----
hdr={}; samples=[]
for line in open(os.path.join(D,"goldbach_sx_checkpoint.csv")):
    if line.startswith("# ") and "=" in line:
        k,v=line[2:].strip().split("=",1); hdr[k]=v
    elif line[0].isdigit():
        N,p,q,src=line.strip().split(","); samples.append((int(N),int(p),int(q),src))
lo,hi=int(hdr["original_start"]),int(hdr["original_end"])
expect=(hi-lo)//2+1
print("CHECKPOINT")
print(f"  version={hdr['cfg_version']}  range=[{lo},{hi}]  last_verified={hdr['last_verified']}")
print(f"  total_verified={hdr['total_verified']}  expected evens={expect}  match={int(hdr['total_verified'])==expect==int(hdr['total_total'])}")
print(f"  total_misses={hdr['total_misses']} unresolved={hdr['unresolved']}")
qh=int(hdr['total_qhot']); cs=int(hdr['cold_survivors']); mr=int(hdr['mr_calls'])
print(f"  ring hit rate={100*qh/expect:.6f}%  mr_calls==cold_survivors: {mr==cs}")
bad=[]
for N,p,q,src in samples:
    ok = p+q==N and isprime(p) and isprime(p+6) and isprime(q) and lo<=N<=hi and p<=int(hdr['cfg_p_anchor_limit'])
    if not ok: bad.append((N,p,q,src))
print(f"  samples={len(samples)}  all witnesses valid (p,p+6,q prime, p+q=N): {not bad}  bad={bad[:3]}")
# ---- misses ----
mrows=[l for l in open(os.path.join(D,"goldbach_sx_misses.csv")) if l[0].isdigit()]
print(f"MISSES file: {len(mrows)} rows")
# ---- minp ----
rows=[]
for line in open(os.path.join(D,"goldbach_sx_minp.csv")):
    if line[0].isdigit():
        N,p,q,t,pt,pc=line.strip().split(","); rows.append((int(N),int(p),int(q),t,int(pt),int(pc)))
print(f"MINP rows={len(rows)}")
types=set(r[3] for r in rows); print("  types:",types, " p_twin==p:",all(r[4]==r[1] for r in rows)," p_cousin all 0:",all(r[5]==0 for r in rows))
# validate every row cheaply: p+q==N, p sexy starter, q prime (sympy BPSW)
badm=0
for N,p,q,*_ in rows:
    if p+q!=N or not (isprime(p) and isprime(p+6) and isprime(q)): badm+=1
print(f"  rows with invalid witness: {badm}")
# minimality check on a random subset: no smaller sexy starter p' with N-p' prime
import random; random.seed(1)
sub=random.sample(rows,2000)
sexy=[p for p in range(3,20000) if isprime(p) and isprime(p+6)]
nonmin=0
for N,p,*_ in sub:
    for s in sexy:
        if s>=p: break
        if isprime(N-s): nonmin+=1; break
print(f"  minimality spot-check (2000 rows): non-minimal = {nonmin}")
# per class stats
print("  fuel by N mod 3:")
for c in (0,1,2):
    ps=[r[1] for r in rows if r[0]%3==c]
    print(f"    {c}: n={len(ps):,} mean {statistics.fmean(ps):7.1f}  median {statistics.median(ps):5.0f}  max {max(ps):5d}")
ps=[r[1] for r in rows]; mx=max(rows,key=lambda r:r[1])
print(f"  overall mean {statistics.fmean(ps):.1f}  max {mx[1]} at N={mx[0]:,}")
# class-split predicted by mod 3: for N=1 mod 3 only p=2 mod 3 usable, N=2 only p=1 mod3, N=0 both
for c in (0,1,2):
    rs=[r for r in rows if r[0]%3==c]
    k=[sum(1 for r in rs if r[1]%3==m) for m in (0,1,2)]
    print(f"    N mod 3={c}: minimal p by p mod 3 (0,1,2) = {k}")
