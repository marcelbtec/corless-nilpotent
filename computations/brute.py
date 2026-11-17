# Brute force: sign sequences s in {+1,-1}^(m-1) with p_m = lambda^m, where
# p_0 = 1, p_1 = lambda, p_k = lambda p_{k-1} + s_{k-1} p_{k-2}  (s = u^2, u in {1, i}).
# Polynomials are stored by the coefficients of lambda^(k-2j), j = 0..k//2 (matching sums).
import itertools, sys

def charpoly(s):
    # returns list c[j] = coefficient of lambda^(m-2j), m = len(s)+1
    P0, P1 = [1], [1]          # p_0, p_1 as matching-sum lists
    for k in range(2, len(s) + 2):
        # p_k = lambda p_{k-1} + s_{k-1} p_{k-2}:  c_k[j] = c_{k-1}[j] + s * c_{k-2}[j-1]
        n = k // 2 + 1
        P2 = [(P1[j] if j < len(P1) else 0) + (s[k - 2] * P0[j - 1] if 1 <= j <= len(P0) else 0) for j in range(n)]
        P0, P1 = P1, P2
    return P1

def nilpotent(s):
    return all(c == 0 for c in charpoly(s)[1:])

def family(k):
    # Corless's recursive family at dimension 2^k - 1, as sign sequences (u = 1 -> +1, u = i -> -1)
    fam = [()]
    for _ in range(k - 1):
        fam = [f + (e, -e) + tuple(reversed(f)) for f in fam for e in (1, -1)]
    return set(fam)

if __name__ == '__main__':
  for m in [2, 3, 4, 5, 6, 7, 8, 9, 15]:
      sols = {s for s in itertools.product((1, -1), repeat=m - 1) if nilpotent(s)}
      k = (m + 1).bit_length() - 1
      fam = family(k) if m == 2 ** k - 1 else set()
      print(f"m={m:2d}: {len(sols):3d} nilpotent sequences; Corless family {len(fam):3d}; equal: {sols == fam}")
      if m in (3, 7):
          for s in sorted(sols): print('      ', ''.join('+' if x > 0 else '-' for x in s))
