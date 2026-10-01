#!/usr/bin/env python3
"""Compare the witness samples of two runs over the same range.

Usage: python3 compare_witnesses.py <reference.csv> <rerun.csv>

COLD rows are minimal partitions at block starts, so with the same block_bits
they must be identical between runs. QHOT rows may legitimately differ (a
different valid partition); they are reported for information only.
Only rows inside the rerun's verified range are compared.
Exit status 0 only if the COLD witnesses are identical.
"""
import sys

def load(f):
    h = {l[2:].split('=')[0]: l.split('=')[1].strip() for l in open(f) if l.startswith('# ') and '=' in l}
    rows = [l.strip().split(',') for l in open(f) if l.count(',') == 3 and l[0].isdigit()]
    return h, rows

def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    ha, ra = load(sys.argv[1]); hb, rb = load(sys.argv[2])
    lo, hi = int(hb['original_start']), int(hb['last_verified'])
    if ha.get('cfg_block_bits') != hb.get('cfg_block_bits'):
        print("warning: block_bits differ; COLD rows are not expected to coincide")
    pick = lambda rows, tag: {r[0]: r for r in rows if r[3] == tag and lo <= int(r[0]) <= hi}
    ca, cb, qa, qb = pick(ra, 'COLD'), pick(rb, 'COLD'), pick(ra, 'QHOT'), pick(rb, 'QHOT')
    same_cold = sum(ca.get(k) == v for k, v in cb.items())
    print(f"compared range [{lo}, {hi}]")
    print(f"COLD: reference {len(ca)} | rerun {len(cb)} | identical {same_cold} | same N set {set(ca) == set(cb)}")
    print(f"QHOT: reference {len(qa)} | rerun {len(qb)} | same N {len(set(qa) & set(qb))} | "
          f"identical partition {sum(qa.get(k) == v for k, v in qb.items())} (information only)")
    ok = set(ca) == set(cb) and same_cold == len(cb)
    print("RESULT:", "COLD witnesses identical" if ok else "COLD WITNESSES DIFFER")
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
