"""Parallel independent check of a DUMP_ALL certificate file (streams; low memory).
Usage: python3 verify_dump_parallel.py <dump.csv> <start> <end> <p_limit> [workers]
Checks every row: N even, p + q == N, p an odd prime <= p_limit (sieve), q prime
(SymPy isprime, deterministic below 2^64), and that every even N in [start, end]
appears exactly once, in order."""
import sys, os, collections
from multiprocessing import Pool
from sympy import isprime

def check_q(chunk):
    return [q for q in chunk if not isprime(q)]

def main():
    path, start, end, plim = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else os.cpu_count()
    sieve = bytearray([1]) * (plim + 1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(plim ** 0.5) + 1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    src = collections.Counter(); bad = []; expect = start; rows = 0
    CH = 20000; batch = []; bad_q = []
    def batches():
        nonlocal expect, rows, batch
        with open(path) as f:
            for line in f:
                N, p, q, s = line.rstrip('\n').split(',')
                N, p, q = int(N), int(p), int(q)
                rows += 1; src[s] += 1
                if N != expect:
                    if len(bad) < 10: bad.append(('gap_or_dup', N, expect))
                    expect = N
                expect += 2
                if N & 1 or p + q != N or not (p & 1) or p > plim or not sieve[p]:
                    if len(bad) < 10: bad.append(('row', N, p, q, s))
                batch.append(q)
                if len(batch) >= CH:
                    yield batch; batch = []
        if batch: yield batch
    with Pool(workers) as pool:
        for res in pool.imap(check_q, batches(), chunksize=4):
            bad_q.extend(res[:10 - len(bad_q)] if len(bad_q) < 10 else [])
            if res and len(bad) < 10: bad.append(('q_not_prime', res[0]))
    if expect != end + 2: bad.append(('range_end', expect - 2, end))
    nfail = len(bad)
    print(f"rows={rows} sources={dict(src)} failures={nfail}")
    for b in bad: print("  ", b)

if __name__ == '__main__':
    main()
