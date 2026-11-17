# Test: if D(s) = 1 mod 2^r at dimension m = 2N+1, how well do the halves do?
# level(X) = largest r' <= CAP with D(X) = 1 mod 2^r' (CAP means exactly 1).
from brute import charpoly
from levels import params
CAP = 10
def level(s):
    c = charpoly(list(s))
    r = 0
    while r < CAP and all(v % (1 << (r + 1)) == 0 for v in c[1:]): r += 1
    return r
for K in (4, 5):
    m, n, t0, basis = params(K)
    N = (m - 1) // 2
    stats = {}
    for i in range(1 << n):
        t = t0
        for b in range(n):
            if (i >> b) & 1: t ^= basis[b]
        s = tuple(-1 if (t >> e) & 1 else 1 for e in range(m - 1))
        r = level(s)
        sig, tau = s[:N - 1], s[N + 1:]           # left half: vertices 1..N, right half: vertices N+2..m
        mid = (s[N - 1], s[N])
        key = r
        st = stats.setdefault(key, {'count': 0, 'minL': CAP, 'minR': CAP, 'mid_opposite': 0, 'tau_rev': 0})
        st['count'] += 1
        st['minL'] = min(st['minL'], level(sig)); st['minR'] = min(st['minR'], level(tau))
        st['mid_opposite'] += mid[0] == -mid[1]
        st['tau_rev'] += tau == tuple(reversed(sig))
    print(f'--- m={m} (halves have dimension {N})')
    for r in sorted(stats):
        st = stats[r]
        print(f'  exact level r={r if r < CAP else "inf"}: {st["count"]:6d} seqs | min level of left half {st["minL"]}, right half {st["minR"]} | middle opposite {st["mid_opposite"]} | right = reversed left {st["tau_rev"]}')
