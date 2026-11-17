# Figure: the 16 nilpotent vectors u in {1, i}^30 at m = 31, as a grid with the recursion's levels marked.
# Usage: python3 make_figure.py   (writes figures/family_31.png)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from brute import family

SURFACE, INK, INK2 = '#fcfcfb', '#0b0b0b', '#52514e'
COL_ONE, COL_I = '#2a78d6', '#eb6834'   # u_k = 1 (sign +1), u_k = i (sign -1)

rows = sorted(family(5))               # sign sequences of length 30, +1 for u = 1, -1 for u = i
assert len(rows) == 16 and all(len(r) == 30 for r in rows)

PAIR_GAP = 0.45                        # extra space between the 15 pairs (2j, 2j+1)
def x0(c): return c + (c // 2) * PAIR_GAP
W = x0(29) + 1

fig, ax = plt.subplots(figsize=(9.6, 5.8), dpi=200)
fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
g = 0.06
for r, s in enumerate(rows):
    for c, v in enumerate(s):
        ax.add_patch(Rectangle((x0(c) + g, r + g), 1 - 2 * g, 1 - 2 * g,
                               facecolor=COL_ONE if v == 1 else COL_I, edgecolor='none'))
# ruler: the centre of each level of the recursion; every pair (2j, 2j+1) is the middle pair of one level
levels = [([14], 1.3, 1.4), ([6, 22], 0.9, 1.0), ([2, 10, 18, 26], 0.6, 0.8),
          ([0, 4, 8, 12, 16, 20, 24, 28], 0.35, 0.7)]
for cells, h, lw in levels:
    for c in cells:
        x = x0(c) + 1                  # boundary between cell c and c + 1
        ax.plot([x, x], [-0.25 - h, -0.25], color=INK if h > 1 else INK2, lw=lw)
ax.plot([x0(14) + 1] * 2, [-0.25, 16.1], color=INK, lw=1.0)
ax.set_xlim(-0.3, W + 0.3); ax.set_ylim(16.9, -3.2)
ax.set_aspect('equal'); ax.axis('off')
ax.text(0, -2.4, 'The 16 nilpotent matrices of dimension 31: superdiagonals u₁ … u₃₀, one per row',
        ha='left', va='bottom', fontsize=10.5, color=INK)
ax.legend(handles=[Patch(color=COL_ONE, label='uₖ = 1'), Patch(color=COL_I, label='uₖ = i')],
          loc='upper left', bbox_to_anchor=(0, 0.02), ncol=2, frameon=False, fontsize=9,
          labelcolor=INK, handlelength=1.2)
fig.tight_layout()
fig.savefig('figures/family_31.png', facecolor=SURFACE)
print('wrote figures/family_31.png')
