# 2-adic reduction of nilpotency. With s_e = 1 - 2 t_e (t_e in {0,1}), the coefficient of x^j in the
# matching polynomial D(x) = sum_M prod_{e in M} s_e x^|M| of the path P_m is
#   sum_S (-2)^|S| prod_{e in S} t_e * N_j(S),   N_j(S) = # j-matchings containing the matching S.
# Mod 2 this is C(m-j, j) (even for all j>=1 iff m = 2^K - 1). Mod 4 it gives LINEAR equations over F_2:
#   sum_e N_j({e}) t_e = C(m-j, j)/2  (mod 2),  j = 1..n.
# Prints the rank of that system and hence the size of the remaining search space.
from math import comb
from brute import family

def nmatch(v, i):          # i-matchings in a path with v vertices
    return comb(v - i, i) if 0 <= i and 2 * i <= v else (1 if (v <= 0 and i == 0) else 0)

def N_edge(m, j, e):        # j-matchings of P_m containing edge e (between vertices e and e+1), e = 1..m-1
    left, right = e - 1, m - e - 1
    return sum(nmatch(left, i) * nmatch(right, j - 1 - i) for i in range(j))

def system(m):
    n = (m - 1) // 2
    rows = []
    for j in range(1, n + 1):
        row = 0
        for e in range(1, m):
            if N_edge(m, j, e) & 1: row |= 1 << (e - 1)
        rhs = (comb(m - j, j) // 2) & 1
        rows.append((row, rhs))
    return rows

def rank_and_solve(rows, nvar):
    piv = {}
    for row, rhs in rows:
        for c in sorted(piv, reverse=True):
            pass
        r, b = row, rhs
        for c in range(nvar - 1, -1, -1):
            if not (r >> c) & 1: continue
            if c in piv: r ^= piv[c][0]; b ^= piv[c][1]
            else: piv[c] = (r, b); break
        else:
            if b: return None, None   # inconsistent
    return len(piv), piv

if __name__ == '__main__':
    for K in range(2, 9):
        m = 2 ** K - 1
        rows = system(m)
        rk, piv = rank_and_solve(rows, m - 1)
        free = (m - 1) - rk
        # check the family satisfies the system
        ok = True
        if K <= 6:
            for s in family(K):
                t = sum(((1 - v) // 2) << i for i, v in enumerate(s))
                ok &= all((bin(row & t).count('1') & 1) == rhs for row, rhs in rows)
        print(f'm={m:3d}: unknowns {m-1:3d}, equations {len(rows):3d}, rank {rk:3d}, free {free:3d} -> {2**free:.3g} candidates; family satisfies: {ok}')
