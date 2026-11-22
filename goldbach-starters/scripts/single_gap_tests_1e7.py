import numpy as np, time
L=10_000_000
s=np.ones(L+20,bool); s[:2]=False
for i in range(2,int((L+20)**.5)+1):
    if s[i]: s[i*i::i]=False
def run(gaps):
    starters=[r for r in np.nonzero(s[:L])[0] if r>=3 and any(s[r+g] for g in gaps)]
    unc=np.arange(6,L+1,2,dtype=np.int64)
    for r in starters:
        m=unc>r
        if not m.any(): break
        hit=np.zeros(len(unc),bool); hit[m]=s[unc[m]-r]
        unc=unc[~hit]
        if len(unc)==0: break
    return unc
for gaps,name in [((2,),"twin only"),((4,),"cousin only"),((2,4),"twin or cousin"),((6,),"sexy only (p+6)"),((8,),"gap 8 only"),((12,),"gap 12 only")]:
    t=time.time(); u=run(gaps)
    by=[int(np.count_nonzero(u%3==c)) for c in (0,1,2)]
    print(f"{name:18} failures to 1e7: {len(u):>9,}  by N mod 3 (0,1,2)={by}  first: {u[:6].tolist()}  ({time.time()-t:.0f}s)")
