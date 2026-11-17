# Affine parametrisation of the mod-4 constraints: t = t0 xor (span of basis), as 64-bit masks
# (bit e-1 = t_e, t_e = 1 means s_e = -1, i.e. u_e = i). Usage: python3 gen_params.py K > params_K.txt
import sys
from mod4 import system
K = int(sys.argv[1]); m = 2 ** K - 1; nvar = m - 1
rows = system(m)
# Gaussian elimination to reduced row echelon form over F_2
piv_rows = []   # list of (pivot col, row mask, rhs)
for row, rhs in rows:
    r, b = row, rhs
    for (c, pr, pb) in piv_rows:
        if (r >> c) & 1: r ^= pr; b ^= pb
    if r == 0:
        assert b == 0, 'inconsistent'
        continue
    c = r.bit_length() - 1
    # eliminate c from existing rows
    piv_rows = [((pc, pr ^ r, pb ^ b) if (pr >> c) & 1 else (pc, pr, pb)) for (pc, pr, pb) in piv_rows]
    piv_rows.append((c, r, b))
pivots = {c for c, _, _ in piv_rows}
free = [c for c in range(nvar) if c not in pivots]
t0 = 0
for c, r, b in piv_rows:
    if b: t0 |= 1 << c          # free vars = 0
basis = []
for f in free:
    v = 1 << f
    for c, r, b in piv_rows:
        if (r >> f) & 1: v |= 1 << c
    basis.append(v)
# sanity: every basis vector solves the homogeneous system, t0 the inhomogeneous one
for row, rhs in rows:
    assert bin(row & t0).count('1') % 2 == rhs
    for v in basis: assert bin(row & v).count('1') % 2 == 0
print(m, len(basis)); print(f'{t0:x}')
for v in basis: print(f'{v:x}')
