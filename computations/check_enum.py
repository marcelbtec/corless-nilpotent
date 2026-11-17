# Compares the output of the exhaustive search (sol_K.txt, one hexadecimal mask per line, bit e-1 set
# when s_e = -1, i.e. u_e = i) with Corless's recursive family at dimension m = 2^K - 1.
# Usage: python3 check_enum.py   (checks every sol_K.txt present, K = 3..6)
import os, sys
from brute import family, nilpotent

ok = True

for K in range(3, 7):
    path = f'sol_{K}.txt'
    if not os.path.exists(path):
        continue
    m = 2 ** K - 1
    masks = [int(w, 16) for w in open(path).read().split()]
    found = {tuple(-1 if (t >> i) & 1 else 1 for i in range(m - 1)) for t in masks}
    fam = family(K)
    ok &= found == fam and all(nilpotent(s) for s in found)
    print(f'm={m:2d}: {len(masks):2d} solutions, all nilpotent: {all(nilpotent(s) for s in found)}; '
          f'equal to the family ({len(fam)}): {found == fam}')
sys.exit(0 if ok else 1)
