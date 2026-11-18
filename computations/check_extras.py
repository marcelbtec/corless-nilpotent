# Exact checks of the statements of the paper that are not formalised in Lean:
#   the characteristic polynomial of J_m(a, c) and the diagonal similarity to A_m^b, the
#   classification of nilpotent J_m(a, c) and tridiagonal sign matrices, the 2^N solutions modulo 4,
#   the single Jordan block, the continued fraction and the resolvent of a nilpotent vector,
#   the doubling identity and its bijection, the sequence c~, the population {0, 1, i}, and the
#   multiplicities of the eigenvalue 0 at m = 15.
# Characteristic polynomials are computed independently of the three-term recurrence, by exact
# Gaussian elimination over the Gaussian rationals at m + 1 points.
# Usage: python3 check_extras.py      (standard library only; a few seconds)
import itertools, random, sys
from fractions import Fraction
from brute import charpoly, family

random.seed(20261117)
ok_all = True
def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f'{name}: {ok}')

# ---- Gaussian rationals as pairs (re, im) of Fractions -------------------------------------
def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
def gsub(a, b): return (a[0] - b[0], a[1] - b[1])
def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def ginv(a):
    n = a[0] * a[0] + a[1] * a[1]
    return (a[0] / n, -a[1] / n)
def gdiv(a, b): return gmul(a, ginv(b))
ZERO, ONE, I = (Fraction(0), Fraction(0)), (Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))
def gq(x): return (Fraction(x.real if isinstance(x, complex) else x), Fraction(x.imag if isinstance(x, complex) else 0))

def det(M):
    M = [row[:] for row in M]; n = len(M); d = ONE
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != ZERO), None)
        if p is None: return ZERO
        if p != c: M[c], M[p] = M[p], M[c]; d = gsub(ZERO, d)
        d = gmul(d, M[c][c])
        for r in range(c + 1, n):
            if M[r][c] != ZERO:
                f = gdiv(M[r][c], M[c][c])
                M[r] = [gsub(M[r][k], gmul(f, M[c][k])) for k in range(n)]
    return d

def J(m, a, c):
    M = [[ZERO] * m for _ in range(m)]
    for k in range(m - 1): M[k][k + 1] = a[k]; M[k + 1][k] = c[k]
    return M

def charpoly_matches(M, s):
    """det(lam I - M) == lam^m D_s(lam^-2) at m + 1 integer points (a polynomial identity of degree m)."""
    m = len(M); d = charpoly(list(s))        # d[j] = coefficient of x^j in D_s
    for lam in range(1, m + 2):
        L = [[gsub((Fraction(lam) if i == j else Fraction(0), Fraction(0)), M[i][j]) for j in range(m)] for i in range(m)]
        rhs = sum(d[j] * Fraction(lam) ** (m - 2 * j) for j in range(len(d)))
        if det(L) != (rhs, Fraction(0)): return False
    return True

def rand_gauss():
    while True:
        z = (Fraction(random.randint(-3, 3)), Fraction(random.randint(-3, 3)))
        if z != ZERO: return z

# ---- Characteristic polynomial of J_m(a, c) with a_k c_k = +-1, and the diagonal similarity -----
ok = True
for _ in range(60):
    m = random.randint(1, 9)
    a = [rand_gauss() for _ in range(m - 1)]
    prods = [random.choice((1, -1)) for _ in range(m - 1)]
    c = [gdiv((Fraction(p), Fraction(0)), a[k]) for k, p in enumerate(prods)]
    ok &= charpoly_matches(J(m, a, c), [-p for p in prods])
report('charpoly of J_m(a,c): det(lam I - J) = lam^m D_s(lam^-2), s = -ac (60 random a, c)', ok)

ok = True
for _ in range(60):
    m = random.randint(2, 10)
    a = [rand_gauss() for _ in range(m - 1)]
    prods = [random.choice((1, -1)) for _ in range(m - 1)]
    c = [gdiv((Fraction(p), Fraction(0)), a[k]) for k, p in enumerate(prods)]
    delta = [ONE]
    for k in range(m - 1): delta.append(gdiv(delta[k], a[k]))
    A = J(m, [ONE] * (m - 1), [(Fraction(p), Fraction(0)) for p in prods])
    S = [[gdiv(gmul(delta[i], A[i][j]), delta[j]) for j in range(m)] for i in range(m)]
    ok &= S == J(m, a, c)
report('diagonal similarity: Delta A_m^b Delta^-1 = J_m(a, c), b = ac (60 random a, c)', ok)

# ---- J_m(a, c) is nilpotent iff m = 2^K - 1 and (a_k c_k) in F_K, m <= 9 -----------------------
def is_nilpotent_charpoly(M):
    m = len(M)
    for lam in range(1, m + 2):
        L = [[gsub((Fraction(lam) if i == j else Fraction(0), Fraction(0)), M[i][j]) for j in range(m)] for i in range(m)]
        if det(L) != (Fraction(lam) ** m, Fraction(0)): return False
    return True
fams = {K: family(K) for K in range(1, 9)}
ok = True
for m in range(1, 10):
    K = (m + 1).bit_length() - 1
    target = fams[K] if m == 2 ** K - 1 else set()
    for prods in itertools.product((1, -1), repeat=m - 1):
        a = [rand_gauss() for _ in range(m - 1)]
        c = [gdiv((Fraction(p), Fraction(0)), a[k]) for k, p in enumerate(prods)]
        ok &= is_nilpotent_charpoly(J(m, a, c)) == (tuple(prods) in target)
report('J_m(a,c) nilpotent iff (a_k c_k) in F_K, all sign patterns, m <= 9 (random a)', ok)
report('F_K closed under negation, K <= 8', all({tuple(-v for v in s) for s in fams[K]} == fams[K] for K in fams))

# ---- D_s = 1 mod 4 has exactly 2^N solutions, all folded, m = 7 and 15 ------------------------
ok = True
for K in (3, 4):
    m = 2 ** K - 1; N = (m - 1) // 2
    sols = [s for s in itertools.product((1, -1), repeat=m - 1)
            if all(c % 4 == (1 if j == 0 else 0) for j, c in enumerate(charpoly(list(s))))]
    folded = all(s[N - 1] == -s[N] and all(s[m - e - 1] == s[e - 1] for e in range(1, N)) for s in sols)
    ok &= len(sols) == 2 ** N and folded
report('D_s = 1 mod 4: exactly 2^N solutions, all folded (m = 7, 15)', ok)

# ---- Continued fraction, resolvent and Jordan block: (T^k)_11 by exact row iteration, K <= 8 ---
def first_row_powers(u):
    """Yield (k, (T^k)_{1,*}) for k = 0..m, T = T_m(u) with Gaussian-integer entries (complex ints)."""
    m = len(u) + 1
    row = [0] * m; row[0] = 1
    for k in range(m + 1):
        yield k, row
        new = [0] * m
        for j in range(m):
            if row[j] == 0: continue
            if j + 1 < m: new[j + 1] += row[j] * u[j]       # T[j][j+1] = u_j
            if j >= 1: new[j - 1] += row[j] * (-u[j - 1])   # T[j][j-1] = -u_{j-1}
        row = new

ok_cf = ok_res = ok_jordan = ok_prod = True
for K in range(2, 9):
    m = 2 ** K - 1
    allowed = {2 ** j - 2 for j in range(1, K + 1)}
    for s in fams[K]:
        eps = {j: s[2 ** (j - 1) - 2] for j in range(2, K + 1)}         # eps_j = s_{2^(j-1)-1}
        want = [0] * (2 ** (K - 1)); want[0] = 1; want[1] = -eps[2]
        for j in range(3, K + 1): want[2 ** (j - 1) - 1] = eps[j]
        got = charpoly(list(s[1:]))
        got = got + [0] * (len(want) - len(got))
        ok_cf &= got[:len(want)] == want and all(v == 0 for v in got[len(want):])
        u = [1 if v == 1 else 1j for v in s]
        coeff = {0: 1, 2: -eps[2]}
        coeff.update({2 ** j - 2: eps[j] for j in range(3, K + 1)})
        for k, row in first_row_powers(u):
            if k < m: ok_res &= row[0] == coeff.get(k, 0)
            if k == m - 1: ok_jordan &= any(v != 0 for v in row)     # T^(m-1) != 0
            if k == m: ok_jordan &= all(v == 0 for v in row)        # first row of T^m is 0
        if K >= 3:
            sig = s[:2 ** (K - 1) - 2]
            p = 1
            for v in sig: p *= v
            ok_prod &= p == -1
report('continued fraction: D_{s\'} = 1 - eps_2 x + sum eps_j x^(2^(j-1)-1), all s in F_K, K <= 8', ok_cf)
report('resolvent: (T^k)_11 = its coefficients, zero off k = 2^j - 2, K <= 8', ok_res)
report('single Jordan block: T^(m-1) != 0, first row of T^m zero, K <= 8', ok_jordan)
report('Product of the signs of sigma in F_{K-1} is -1 (3 <= K <= 8)', ok_prod)

ok = True
for K in range(2, 9):
    N = 2 ** (K - 1) - 1
    for s in fams[K]:
        sig, eps = s[:N - 1], s[N - 1]
        p = 1
        for v in sig: p *= v
        lhs = charpoly(list(s[1:]))
        rhs = charpoly(list(sig[1:])) if N > 1 else [1]
        rhs = rhs + [0] * (N + 1 - len(rhs)); rhs[N] += -eps * p
        trim = lambda a: a[:max(i for i, v in enumerate(a) if v != 0) + 1]
        ok &= trim(lhs) == trim(rhs)
report('recursion step: D_{s\'} = D_{sigma\'} - eps prod(sigma) x^N, all s in F_K, K <= 8', ok)

# ---- Doubling identity: D_(th1,-th1,...,thN,-thN)(x) = D_tau(x^2) -----------------------------
ok = True
for _ in range(500):
    N = random.randint(1, 40)
    th = [random.choice((1, -1)) for _ in range(N)]
    r = [v for t in th for v in (t, -t)]
    tau = [th[j] * th[j + 1] for j in range(N - 1)]
    lhs = charpoly(r); dt = charpoly(tau)
    rhs = [0] * (2 * len(dt) - 1)
    for j, v in enumerate(dt): rhs[2 * j] = v
    trim = lambda a: a[:max([i for i, v in enumerate(a) if v != 0] + [0]) + 1]
    ok &= trim(lhs) == trim(rhs)
report('doubling identity D_r(x) = D_tau(x^2), 500 random theta, N <= 40', ok)

# ---- Squares: (A_{2N+1}^{-r})^2 acts on the even-numbered basis vectors as A_N^{-tau} ----------
ok = True
for _ in range(200):
    N = random.randint(1, 12)
    th = [random.choice((1, -1)) for _ in range(N)]
    r = [v for x in th for v in (x, -x)]
    tau = [th[j] * th[j + 1] for j in range(N - 1)]
    n = 2 * N + 1
    M = [[0] * n for _ in range(n)]
    for k in range(n - 1): M[k][k + 1] = 1; M[k + 1][k] = -r[k]
    M2 = [[sum(M[i][k] * M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    ev = [i for i in range(n) if i % 2 == 1]          # 0-based positions of the vertices 2, 4, ..., 2N
    od = [i for i in range(n) if i % 2 == 0]
    want = [[0] * N for _ in range(N)]
    for k in range(N - 1): want[k][k + 1] = 1; want[k + 1][k] = -tau[k]
    ok &= [[M2[i][j] for j in ev] for i in ev] == want and all(M2[i][j] == 0 for i in od for j in ev)
report('squares: (A_{2N+1}^{-r})^2 acts on the even-numbered basis vectors as A_N^{-tau} (200 random theta)', ok)

# ---- Doubling rule: Phi is a bijection {+-1} x F_{K-1} -> F_K ---------------------------------
def Phi(eta, tau):
    th = [eta]
    for t in tau: th.append(t * th[-1])
    return tuple(v for t in th for v in (t, -t))
ok = all({Phi(eta, tau) for eta in (1, -1) for tau in fams[K - 1]} == fams[K] for K in range(2, 9))
report('doubling bijection Phi({+-1} x F_{K-1}) = F_K, K <= 8', ok)

# ---- The sequence c~ ---------------------------------------------------------------------------
ct = {1: 1}
for j in range(1, 2 ** 9):
    ct[2 * j] = ct[2 * j - 1] * ct[j]
    ct[2 * j + 1] = -ct[2 * j]
ok = True
for K in range(1, 9):
    m = 2 ** K - 1
    b = tuple(ct[i] for i in range(2, m + 1))
    ok &= b in fams[K] and all(v == 0 for v in charpoly([-v for v in b])[1:])
    if K >= 2:
        bprev = tuple(ct[i] for i in range(2, 2 ** (K - 1)))
        ok &= tuple(-v for v in b) == Phi(-1, tuple(-v for v in bprev))
report('sequence c~: b^(K) in F_K, A_m^b nilpotent, -b^(K) = Phi_-1(-b^(K-1)), K <= 8', ok)

# ---- Population {0, 1, i}: nilpotent iff every maximal run of nonzero signs lies in some F_K -----
def runs(s):
    out, cur = [], []
    for v in s:
        if v == 0: out.append(tuple(cur)); cur = []
        else: cur.append(v)
    out.append(tuple(cur)); return out
allfam = set().union(*fams.values())
count = 0; ok = True
for m in range(1, 11):
    for s in itertools.product((0, 1, -1), repeat=m - 1):
        count += 1
        ok &= all(c == 0 for c in charpoly(list(s))[1:]) == all(r in allfam for r in runs(s))
report(f'population {{0, 1, i}}: nilpotent iff all runs in F_K, all {count} sign sequences, m <= 10', ok and count == 29524)

# ---- Multiplicity of the eigenvalue 0 at m = 15: m - 2 deg D_s takes every odd value except 9 ----
mults = set()
for s in itertools.product((1, -1), repeat=14):
    d = charpoly(list(s))
    mults.add(15 - 2 * max(j for j, v in enumerate(d) if v != 0))
report('m = 15: multiplicities of the eigenvalue 0 are 1, 3, 5, 7, 11, 13, 15 (9 does not occur)', mults == {1, 3, 5, 7, 11, 13, 15})

print('ALL CHECKS PASSED' if ok_all else 'SOME CHECK FAILED')
sys.exit(0 if ok_all else 1)
