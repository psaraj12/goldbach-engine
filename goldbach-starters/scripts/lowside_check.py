import sys, time
LIMIT = int(sys.argv[1])
KMAX = 4

t0 = time.time()
size = LIMIT + 2*KMAX + 2
sieve = bytearray([1]) * size
sieve[0:2] = b'\x00\x00'
for i in range(2, int(size**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))

nbytes = (LIMIT + 8) // 8 + 1
buf = bytearray(nbytes)
primes = [i for i in range(2, LIMIT + 1) if sieve[i]]
for p in primes:
    buf[p >> 3] |= 1 << (p & 7)
P_int = int.from_bytes(bytes(buf), 'little')
mask_all = (1 << (LIMIT + 1)) - 1

ebuf = bytearray(nbytes)
for N in range(6, LIMIT + 1, 2):
    ebuf[N >> 3] |= 1 << (N & 7)
E = int.from_bytes(bytes(ebuf), 'little')

uncovered = E
for k in (1, 2):
    if uncovered == 0: break
    tk = time.time()
    Sk = [r for r in primes if sieve[r + 2*k]]
    reach = 0
    for i, r in enumerate(Sk):
        # keep only partners p >= r  ->  N = r + p >= 2r  (r is the lower-or-equal member)
        reach |= (P_int >> r) << (2*r)   # bits >= r of P, shifted by r: bit N set iff N-r prime and N-r >= r
        if i % 4000 == 3999:
            reach &= mask_all
    reach &= mask_all
    newly = uncovered & reach
    uncovered &= ~reach
    print(f"LOW-SIDE step +{2*k}: newly={bin(newly).count('1'):,}  remaining={bin(uncovered).count('1'):,}  ({time.time()-tk:.1f}s)")

if uncovered:
    rem = []
    u = uncovered
    while u and len(rem) < 60:
        n = (u & -u).bit_length() - 1
        rem.append(n); u &= u - 1
    cnt = bin(uncovered).count('1')
    print(f"COUNTEREXAMPLES to low-side claim (need upper-side pair-prime): {cnt:,}")
    print("first ones:", rem)
else:
    print(f"CLAIM HOLDS: every even in [6,{LIMIT:,}] has a lower-half partition prime r (r <= N/2) with r+2 or r+4 also prime.")
print(f"total {time.time()-t0:.1f}s")
