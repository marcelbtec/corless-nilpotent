# Numerical checks of every step of the proof.
import random
from math import comb
from brute import charpoly          # matching sums of the path, i.e. coefficients of D_s

def D(s): return charpoly(list(s))   # D_s(x) as list of coefficients of x^0, x^1, ...
def square(a):
    out = [0] * (2 * len(a) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(a): out[i + j] += u * v
    return out
def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a

# (1) Identity D(sigma, e, -e, rev sigma) = D(sigma)^2 for random sigma of many lengths (not only nilpotent ones)
ok = True
for _ in range(3000):
    L = random.randint(0, 40)
    sig = tuple(random.choice((1, -1)) for _ in range(L)); e = random.choice((1, -1))
    s = sig + (e, -e) + tuple(reversed(sig))
    ok &= trim(D(s)) == trim(square(D(sig)))
print('(1) D(sigma, e, -e, rev sigma) == D(sigma)^2 for 3000 random sigma:', ok)

# (2) Over F_2, f_k(y + y^2) = (1+y)^(k+1) + y^(k+1), where f_k = sum_i C(k-i, i) x^i  (polynomials as int bitmasks)
def pmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r
def ppow(a, k):
    r = 1
    for _ in range(k): r = pmul(r, a)
    return r
from functools import lru_cache
@lru_cache(maxsize=None)
def f2(k):   # f_k mod 2 in x, via f_k = f_(k-1) + x f_(k-2), f_0 = f_1 = 1
    if k <= 1: return 1
    return f2(k - 1) ^ (f2(k - 2) << 1)
def subst(px):  # substitute x = y + y^2
    r, xp, z = 0, 1, 0b110
    i = 0
    while px >> i:
        if (px >> i) & 1: r ^= xp
        xp = pmul(xp, z); i += 1
    return r
ok2 = all(subst(f2(k)) == ppow(0b11, k + 1) ^ (1 << (k + 1)) for k in range(0, 200))
print('(2) f_k(y + y^2) == (1+y)^(k+1) + y^(k+1) over F_2 for k < 200:', ok2)

# (3) The columns c_e = x f_{e-1} f_{m-1-e} (mod 2), e = 1..N, are linearly independent over F_2, m = 2^K - 1
def rank(vecs):
    basis = []
    for v in vecs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
import sys
sys.setrecursionlimit(20000)
for k in range(5000): f2(k)
for K in range(2, 12):
    m = 2 ** K - 1; N = (m - 1) // 2
    cols = [pmul(0b10, pmul(f2(e - 1), f2(m - 1 - e))) for e in range(1, N + 1)]
    mirror = True  # c_e and c_(m-e) are the same product f_(e-1) f_(m-1-e) with the factors swapped
    print(f'(3) K={K:2d} m={m:5d}: rank of c_1..c_N = {rank(cols)} (N = {N})')

# (4) f_N = 1 mod 2 for N = 2^(K-1) - 1
print('(4) f_(2^j - 1) == 1 mod 2 for j < 12:', all(f2(2 ** j - 1) == 1 for j in range(1, 12)), '| and f_k != 1 mod 2 for other k < 4096:', all(f2(k) != 1 for k in range(2, 4096) if (k + 1) & k))
