# Computations

These scripts check every step of the proof in `../proof.md` numerically. They also extend Corless's exhaustive verification from m = 31 to m = 63. `run_all.sh` reruns all of them.

All scripts work with signs s_k = u_k², where s_k = +1 for u_k = 1 and s_k = −1 for u_k = i. T_m(u) is nilpotent exactly when the matching polynomial D_s of the path with edge signs s equals 1.

## Scripts

| Script | What it checks | Output |
|---|---|---|
| `brute.py` | brute force over all sign sequences for m ≤ 9 and m = 15; also the helpers `charpoly` and `family` used by the other scripts | the nilpotent sets equal Corless's family (2, 4 and 8 at m = 3, 7, 15), none at other m |
| `check_proof.py` | Lemma 1 (folding) for 3,000 random σ; the closed form of Lemma 2 for k < 200; Lemma 3 (rank N) for m = 2^K − 1 up to 2,047; parity of f_k for k < 4,096 | all hold |
| `check_matrices.py` | the matrices themselves, in exact Gaussian-integer arithmetic | T^m = 0 and T^(m−1) ≠ 0 for every family member at m = 3, 7, 15; no other nilpotent at m = 7 or among 200 random vectors at m = 15; det(λI − T) equals the recurrence |
| `mod4.py` | the reduction modulo 4: the linear system over 𝔽₂ and its rank, for m up to 255 | full rank; 2^31 candidates remain at m = 63 |
| `gen_params.py` | writes the solution space of the modulo-4 system as an affine space, `params_K.txt` | 64-bit masks |
| `enum.c` | exhaustive search of that affine space; each candidate is filtered by evaluating D_s at two fixed points modulo 2^61 − 1 and then checked exactly | `sol_K.txt` (solutions) and `log_K.txt` (counts per slice) |
| `check_enum.py` | compares `sol_K.txt` with Corless's family | equal for m = 7, 15, 31, 63 (4, 8, 16, 32 solutions) |
| `halves.py`, `levels.py` | exploratory: how far D_s = 1 holds modulo higher powers of 2, and how the two halves behave | tables |
| `sfraction.py` | the continued-fraction observation: for every family member, D_{s₂…} has the form 1 ± x ± x³ ± ⋯ | printed polynomials |
| `check_extras.py` | the statements of the paper that are not formalised in Lean: tridiagonal sign matrices J_m(a, c) and the diagonal similarity, the 2^N solutions modulo 4, the single Jordan block, the continued fraction and the resolvent, the doubling identity and the sequence c̃; characteristic polynomials by exact elimination over the Gaussian rationals | all hold (K ≤ 8 where the statement is about F_K) |
| `make_figure.py` | draws `figures/family_31.png`, the 16 nilpotent vectors at m = 31 | PNG; needs matplotlib |

## Data

- **`params_K.txt`**: the modulo-4 solution space at m = 2^K − 1. The first line gives m and the dimension n. The second line gives an offset t₀ in hexadecimal, followed by n basis vectors, one per line. Bit e − 1 of a mask is t_e, and t_e = 1 means s_e = −1.
- **`sol_K.txt`**: the nilpotent sign sequences found by `enum.c`, one mask per line in the same encoding.
- **`log_K.txt`**: the per-slice counts of candidates, of candidates passing the filter, and of nilpotent sequences.

The search is complete without relying on the proof. Every nilpotent sequence satisfies D_s ≡ 1 (mod 4), so it lies in the searched space.

Usage of the enumerator: `./enum params_K.txt slice nslices`. It prints the solutions of one slice to stdout and the counts to stderr.
