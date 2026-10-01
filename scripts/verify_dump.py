"""Independent check of a v5 DUMP_ALL certificate file.
For every row: N even, p + q == N, p an odd prime <= p_limit, q prime (SymPy,
deterministic below 2^64). Also checks every even N in [start, end] appears exactly once."""
import sys, collections
from sympy import isprime
path, start, end, plim = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
sieve = bytearray([1]) * (plim + 1); sieve[0:2] = b'\x00\x00'
for i in range(2, int(plim**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
bad, src, expect, prev = [], collections.Counter(), start, None
for line in open(path):
    N, p, q, s = line.strip().split(','); N, p, q = int(N), int(p), int(q)
    src[s] += 1
    if N != expect: bad.append(('gap_or_dup', N, expect)); expect = N
    expect += 2
    if N % 2 or p + q != N or p % 2 == 0 or p > plim or not sieve[p] or not isprime(q):
        bad.append((N, p, q, s))
if expect != end + 2: bad.append(('range_end', expect - 2, end))
print(f"rows={sum(src.values())} sources={dict(src)} failures={len(bad)}")
for b in bad[:10]: print("  ", b)
