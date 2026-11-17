# Ties the theorem to Corless's matrices: exact Gaussian-integer arithmetic on T_m(u), with
# superdiagonal u (entries 1 or i), subdiagonal -u, zero diagonal. Gaussian integers are pairs (a, b) = a + b i.
import random, itertools
from brute import family, charpoly
def gmul(p, q): return (p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0])
def gadd(p, q): return (p[0] + q[0], p[1] + q[1])
Z = (0, 0)
def T(u):
    m = len(u) + 1
    A = [[Z] * m for _ in range(m)]
    for k, x in enumerate(u): A[k][k + 1] = x; A[k + 1][k] = (-x[0], -x[1])
    return A
def matmul(A, B):
    n = len(A)
    return [[(lambda r, c: __import__('functools').reduce(gadd, (gmul(A[r][k], B[k][c]) for k in range(n)), Z))(r, c) for c in range(n)] for r in range(n)]
def power(A, k):
    R = A
    for _ in range(k - 1): R = matmul(R, A)
    return R
def is_zero(A): return all(v == Z for row in A for v in row)
ONE, I = (1, 0), (0, 1)
def u_of(s): return [ONE if v > 0 else I for v in s]          # s = u^2: u = 1 -> +1, u = i -> -1

# (a) every family member: T^m = 0 and T^(m-1) != 0 (single Jordan block)
for K in (2, 3, 4):
    m = 2 ** K - 1
    res = [(is_zero(power(T(u_of(s)), m)), not is_zero(power(T(u_of(s)), m - 1))) for s in family(K)]
    print(f'(a) m={m}: {len(res)} family members, all T^m = 0: {all(r[0] for r in res)}, all T^(m-1) != 0: {all(r[1] for r in res)}')
# (b) a non-family u at m = 7 and 15 is not nilpotent; exhaustive at m = 7
fam7 = family(3)
bad = [s for s in itertools.product((1, -1), repeat=6) if s not in fam7 and is_zero(power(T(u_of(s)), 7))]
print(f'(b) m=7: non-family vectors with T^7 = 0: {len(bad)} of {64 - len(fam7)}')
random.seed(1); fam15 = family(4)
smp = [tuple(random.choice((1, -1)) for _ in range(14)) for _ in range(200)]
print(f'(b) m=15: 200 random non-family vectors with T^15 = 0: {sum(1 for s in smp if s not in fam15 and is_zero(power(T(u_of(s)), 15)))}')
# (c) characteristic polynomial from the matrix (Faddeev-LeVerrier, exact) equals the recurrence
def charpoly_matrix(A):
    n = len(A); M = [[Z] * n for _ in range(n)]; coeffs = [ONE]; Id = [[ONE if r == c else Z for c in range(n)] for r in range(n)]
    for k in range(1, n + 1):
        M = [[gadd(M2, (coeffs[-1][0] * Id[r][c][0], coeffs[-1][0] * Id[r][c][1]) if True else Z) for c, M2 in enumerate(row)] for r, row in enumerate(matmul(A, M) if k > 1 else [[Z] * n for _ in range(n)])]
        if k == 1: M = Id
        AM = matmul(A, M); tr = Z
        for r in range(n): tr = gadd(tr, AM[r][r])
        assert tr[0] % k == 0 and tr[1] % k == 0
        coeffs.append((-tr[0] // k, -tr[1] // k))
    return coeffs        # coefficients of lambda^n, lambda^(n-1), ..., lambda^0
ok = True
for _ in range(60):
    m = random.randint(2, 12); s = tuple(random.choice((1, -1)) for _ in range(m - 1))
    cm = charpoly_matrix(T(u_of(s)))
    rec = charpoly(list(s))            # coefficient of lambda^(m-2j) is rec[j]
    want = [(0, 0)] * (m + 1)
    for j, v in enumerate(rec): want[2 * j] = (v, 0)
    ok &= cm == want
print('(c) det(lambda I - T) from the matrix equals the recurrence for 60 random u, m <= 12:', ok)
