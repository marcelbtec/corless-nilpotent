# For a nilpotent sign sequence s (length 2n), the S-fraction 1/(1 + s_1 x/(1 + s_2 x/(...)))
# equals the +-1 weighted matching polynomial of the path with vertex 1 removed (the full
# matching polynomial being 1). Print it for the family.
from brute import family, charpoly
def show(c): return ' '.join(f'{v:+d}x^{j}' for j, v in enumerate(c) if v)
for k in (2, 3, 4, 5):
    print(f'--- dimension {2**k - 1}')
    for s in sorted(family(k)):
        F = charpoly(s[1:])   # matching sums of the path on vertices 2..m (2n vertices)
        print(''.join('+' if v > 0 else '-' for v in s), '  F =', show(F))
