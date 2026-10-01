#!/usr/bin/env python3
"""Independently check every witness row in a goldbach checkpoint file.

Usage: python3 verify_witnesses.py <checkpoint.csv> [--minimal]

For every row N,p,q,source: N even, p + q == N, p <= the run's anchor limit,
p and q prime (SymPy isprime, deterministic below 2^64), and N inside the
verified range. With --minimal, COLD rows are also checked to be minimal
partitions (no smaller odd prime p' with N - p' prime).
Exit status 0 only if every check passes.
"""
import sys
from sympy import isprime, primerange

def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path, minimal = sys.argv[1], "--minimal" in sys.argv[2:]
    h = {l[2:].split('=')[0]: l.split('=')[1].strip()
         for l in open(path) if l.startswith('# ') and '=' in l}
    rows = [l.strip().split(',') for l in open(path) if l.count(',') == 3 and l[0].isdigit()]
    lo, hi, P = int(h['original_start']), int(h['last_verified']), int(h['cfg_p_anchor_limit'])
    print(f"range [{lo}, {hi}] | verified {h['total_verified']} of {h['total_total']} | "
          f"misses {h['total_misses']} (unresolved {h['unresolved']}) | anchor limit {P}")
    bad, src, nonmin = [], {}, 0
    small = list(primerange(3, 20000)) if minimal else []
    for N, p, q, s in rows:
        N, p, q = int(N), int(p), int(q)
        src[s] = src.get(s, 0) + 1
        if not (N % 2 == 0 and p + q == N and p <= P and lo <= N <= hi and isprime(p) and isprime(q)):
            bad.append((N, p, q, s))
        if minimal and s == 'COLD':
            for t in small:
                if t >= p:
                    break
                if isprime(N - t):
                    nonmin += 1
                    break
    print(f"rows {len(rows)} {src} | failures {len(bad)}" + (f" | non-minimal COLD {nonmin}" if minimal else ""))
    for b in bad[:10]:
        print("  ", b)
    sys.exit(1 if bad or nonmin else 0)

if __name__ == '__main__':
    main()
