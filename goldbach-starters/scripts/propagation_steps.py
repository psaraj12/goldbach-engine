import sys, time
LIMIT = int(sys.argv[1])
KMAX = 8  # steps +2 .. +16

t0 = time.time()
size = LIMIT + 2*KMAX + 2
sieve = bytearray([1]) * size
sieve[0:2] = b'\x00\x00'
for i in range(2, int(size**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))

# pack primality into a big int (bit n = n is prime), limited to <= LIMIT
nbytes = (LIMIT + 8) // 8 + 1
buf = bytearray(nbytes)
primes = [i for i in range(2, LIMIT + 1) if sieve[i]]
for p in primes:
    buf[p >> 3] |= 1 << (p & 7)
P_int = int.from_bytes(bytes(buf), 'little')
print(f"[{LIMIT:,}] sieve+pack done {time.time()-t0:.1f}s, primes={len(primes):,}")

mask_all = (1 << (LIMIT + 1)) - 1
# even-number mask for N in [6, LIMIT]
ebuf = bytearray(nbytes)
for N in range(6, LIMIT + 1, 2):
    ebuf[N >> 3] |= 1 << (N & 7)
E = int.from_bytes(bytes(ebuf), 'little')
total_evens = bin(E).count('1')

uncovered = E
bands = [(6, 10**4), (10**4, 10**5), (10**5, 10**6), (10**6, 10**7)]
bands = [(lo, hi) for lo, hi in bands if lo < LIMIT]
def band_mask(lo, hi):
    hi = min(hi, LIMIT)
    return ((1 << (hi + 1)) - 1) ^ ((1 << lo) - 1)

for k in range(1, KMAX + 1):
    if uncovered == 0: break
    tk = time.time()
    Sk = [r for r in primes if r + 2*k <= LIMIT + 2*KMAX and sieve[r + 2*k]]
    reach = 0
    cnt = 0
    for r in Sk:
        reach |= P_int << r
        cnt += 1
        if cnt % 4000 == 0:
            reach &= mask_all
    reach &= mask_all
    newly = uncovered & reach
    uncovered &= ~reach
    print(f"step +{2*k}: pair-primes={len(Sk):,}  newly covered={bin(newly).count('1'):,}  "
          f"remaining uncovered={bin(uncovered).count('1'):,}  ({time.time()-tk:.1f}s)")
    for lo, hi in bands:
        bm = band_mask(lo, hi)
        tot = bin(E & bm).count('1')
        got = bin(newly & bm).count('1')
        if tot: print(f"    band [{lo:,},{min(hi,LIMIT):,}]: {got:,}/{tot:,} = {100*got/tot:.3f}%")

if uncovered:
    rem = []
    u = uncovered
    while u and len(rem) < 50:
        n = (u & -u).bit_length() - 1
        rem.append(n); u &= u - 1
    print("STRANDED (no step <= +%d):" % (2*KMAX), rem, "..." if uncovered.bit_length() else "")
else:
    print(f"ALL {total_evens:,} evens in [6, {LIMIT:,}] can step forward. Max minimal step used: see above.")
print(f"total {time.time()-t0:.1f}s")
