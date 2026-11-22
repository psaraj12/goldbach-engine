import numpy as np, time
t0=time.time()
LIMIT = 100_000_000
# sieve
sieve = np.ones(LIMIT+5, dtype=bool); sieve[:2]=False
for i in range(2, int((LIMIT+4)**0.5)+1):
    if sieve[i]: sieve[i*i::i]=False
twins = np.nonzero(sieve[:-2] & sieve[2:])[0]; twins = twins[twins>=3]  # r, r+2 prime
print(f"sieve+twins {time.time()-t0:.0f}s, twins={len(twins):,}")

# Class N=2 mod 3 (even -> N=2 mod 6): can-step-+2 iff N-3 prime (r=3 or r=N-3; any other twin r=2 mod 3 forces 3|q).
# Classes N=0,4 mod 6: verify EVERY one is covered by some twin r (layered exact elimination).
unc = np.concatenate([np.arange(6, LIMIT+1, 6, dtype=np.int64),
                      np.arange(10, LIMIT+1, 6, dtype=np.int64)])
for r in twins:
    m = unc > r
    if not m.any(): break
    hit = np.zeros(len(unc), dtype=bool)
    hit[m] = sieve[unc[m]-r]
    unc = unc[~hit]
    if len(unc)==0: break
print(f"classes 0,1: residual uncovered by any twin = {len(unc)}  ({time.time()-t0:.0f}s)")

bands = [(6,10**4),(10**4,10**5),(10**5,10**6),(10**6,10**7),(10**7,10**8)]
print(f"{'band':>16} {'evens':>12} {'needs +4':>10} {'%':>8}")
for lo,hi in bands:
    N2 = np.arange(lo + ((2-lo)%6), hi+1, 6, dtype=np.int64)   # N=2 mod 6 in (lo,hi]
    N2 = N2[(N2>lo)&(N2<=hi)&(N2>=8)]
    needs = int(np.count_nonzero(~sieve[N2-3]))
    needs += int(np.count_nonzero((unc>lo)&(unc<=hi)))         # residuals, if any
    tot = (hi-lo)//2
    print(f"[{lo:>7,}-{hi:>9,}] {tot:>12,} {needs:>10,} {100*needs/tot:>7.3f}%")
print(f"total {time.time()-t0:.0f}s")
