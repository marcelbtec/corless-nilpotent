# For each dimension m = 2^K - 1: among sign sequences satisfying the mod-4 equations (affine space),
# count those whose path matching polynomial D satisfies D = 1 modulo 2^r, for r = 2, 3, 4, ...
from brute import charpoly, family
def params(K):
    lines = open(f'params_{K}.txt').read().split()
    m, n = int(lines[0]), int(lines[1]); t0 = int(lines[2], 16); basis = [int(v, 16) for v in lines[3:3 + n]]
    return m, n, t0, basis
for K in (3, 4, 5):
    m, n, t0, basis = params(K)
    counts = {}
    for i in range(1 << n):
        t = t0
        for b in range(n):
            if (i >> b) & 1: t ^= basis[b]
        s = [(-1 if (t >> e) & 1 else 1) for e in range(m - 1)]
        c = charpoly(s)
        for r in range(2, 9):
            if all(v % (1 << r) == 0 for v in c[1:]): counts[r] = counts.get(r, 0) + 1
    print(f'm={m:2d}: mod-4 space {1 << n:6d};', ' '.join(f'mod 2^{r}: {counts.get(r, 0)}' for r in range(2, 9)), f'| family {len(family(K))}')
