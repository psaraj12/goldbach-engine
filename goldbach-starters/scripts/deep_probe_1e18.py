import random, time
from sympy import isprime

# small-prime sieve for twin/cousin starters
B = 2_000_000
sieve = bytearray([1]) * B
sieve[0:2] = b'\x00\x00'
for i in range(2, int(B**0.5)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
twins   = [r for r in range(3, B-2) if sieve[r] and sieve[r+2]]   # r, r+2 prime
cousins = [r for r in range(3, B-4) if sieve[r] and sieve[r+4]]   # r, r+4 prime
print(f"starters loaded: {len(twins):,} twins, {len(cousins):,} cousins (r < {B:,})")

random.seed(42)
SAMPLES = 300
for expo in (9, 12, 15, 18):
    base = 10**expo
    t0 = time.time()
    max_r, max_probes, tot_probes = 0, 0, 0
    fails = 0
    step_used = {2:0, 4:0}
    for _ in range(SAMPLES):
        N = base + 2*random.randrange(1, 10**7)   # random even just above 10^expo
        found = None
        probes = 0
        # gear choice by N mod 3 (try the unblocked gear first, then the other)
        gears = [(2, twins), (4, cousins)] if N % 3 != 2 else [(4, cousins), (2, twins)]
        for step, starters in gears:
            for r in starters:
                probes += 1
                if isprime(N - r):          # r <= 2e6 << N/2, firmly low-side
                    found = (step, r); break
            if found: break
        if found is None:
            fails += 1
            print(f"  !! no low-side twin/cousin starter < {B:,} for N={N}")
        else:
            step_used[found[0]] += 1
            max_r = max(max_r, found[1]); max_probes = max(max_probes, probes); tot_probes += probes
    print(f"N ~ 10^{expo}: {SAMPLES} samples | fails={fails} | steps used +2:{step_used[2]} +4:{step_used[4]} | "
          f"max starter r needed={max_r} | avg probes={tot_probes/SAMPLES:.1f} | max probes={max_probes} | {time.time()-t0:.1f}s")
