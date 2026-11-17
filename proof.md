# Nilpotent skew-symmetric tridiagonal Bohemian matrices: a proof of Corless's conjecture

Status: 17 November 2026. The proof is formalised in Lean 4 with Mathlib and checked by machine (folder `lean/`, see "Formal verification" below). The formal version covers the statement about Corless's complex matrices in his own terms and the count of nilpotent matrices. It has not yet been read by anyone other than its authors.

Openness, checked before and after the proof (17 November 2026): the conjecture is stated as open in Corless's paper (2021), his talk of December 2021 ("I thought I had a proof that this construction gave all the nilpotents, but I now think I don't") and his Leuven talk of October 2023 ("experimentally checked to m = 31"). The paper has one citing work (Bohemian Matrix Geometry, 2022), which does not address it. No arXiv paper by Corless since 2023, no arXiv abstract combining nilpotent, tridiagonal and Bohemian, and no MathOverflow or Math Stack Exchange thread was found. A second check on the same day found nothing either. It covered Semantic Scholar, OpenCitations and OpenAlex citations, his talk slides from 2024 to 2026, his CV (February 2026) and Maple Transactions from 2024 to 2026. His slides of 11 October 2023 (KU Leuven, slide 21) still call it a conjecture. The argument below is short, and its author expected it to be easy; it may exist unpublished. A message to Corless is the definitive check.

## The conjecture

For a vector u = (u₁, …, u_{m−1}) with entries in {1, i}, let T_m(u) be the m × m skew-symmetric tridiagonal matrix with zero diagonal, superdiagonal u and subdiagonal −u. Corless (2021) proved that T_m(u) can be nilpotent only when m = 2^K − 1, and that the vectors built recursively by

u ↦ (u, 1, i, rev u)   and   u ↦ (u, i, 1, rev u),

starting from the empty vector at m = 1, give nilpotent matrices at every m = 2^K − 1. His Conjecture 5.1 states that these are the only nilpotent matrices in the family. He verified it exhaustively up to m = 31 and noted that m = 63 would take about 73,000 years by the same method; in 2023 he wrote that it "ought to be easily provable, but I failed on my first try".

Reference: R. M. Corless, Skew-symmetric tridiagonal Bohemian matrices, Maple Transactions 1(2) (2021), Article 14360, doi:10.5206/mt.v1i2.14360.

## Theorem

T_m(u) is nilpotent if and only if m = 2^K − 1 for some K ≥ 1 and u belongs to Corless's recursive family. In particular there are exactly 2^{K−1} nilpotent matrices of dimension 2^K − 1 in the family, and none in other dimensions.

## Notation

Put s_k = u_k² ∈ {+1, −1} (s_k = +1 for u_k = 1, s_k = −1 for u_k = i). Expanding det(λI − T_m) along the last row gives the three-term recurrence

p₀ = 1,   p₁ = λ,   p_k = λ·p_{k−1} + s_{k−1}·p_{k−2},

because the two off-diagonal entries of λI − T next to the diagonal are −u and +u, whose product is −u² = −s. Unrolling the recurrence, p_m is a sum over the **matchings** M of the path graph P_m with vertices 1, …, m and edges e = {e, e+1}:

p_m(λ) = Σ_M (∏_{e∈M} s_e) · λ^{m − 2|M|} = λ^m · D_s(λ^{−2}),   D_s(x) = Σ_M (∏_{e∈M} s_e) · x^{|M|} ∈ ℤ[x].

So T_m(u) is nilpotent if and only if D_s = 1. Corless's recursion, in terms of s, is s ↦ (s, ε, −ε, rev s) with ε = ±1. Write F_1 = {()} and F_K = {(σ, ε, −ε, rev σ) : σ ∈ F_{K−1}, ε = ±1}; then |F_K| = 2^{K−1}, since σ and ε can be read off from the result.

For the unweighted path write f_k(x) = Σ_j C(k − j, j)·x^j, the matching polynomial of P_k (a Fibonacci polynomial), with f₀ = f₁ = 1 and f_k = f_{k−1} + x·f_{k−2}.

## Lemma 1 (folding identity)

For every σ ∈ {±1}^{N−1} and ε = ±1, the sequence s = (σ, ε, −ε, rev σ) of length 2N satisfies D_s = D_σ².

Proof. The path has m = 2N + 1 vertices; its middle vertex c = N + 1 is joined to N by an edge of weight ε and to N + 2 by an edge of weight −ε. Let L be the subpath 1, …, N (weights σ) and R the subpath N + 2, …, 2N + 1 (weights rev σ). The reflection v ↦ 2N + 2 − v maps R onto L and preserves edge weights, so D_R = D_L = D_σ, and R without its first vertex has the same matching polynomial as L without its last vertex, call it D_{L′}. Splitting the matchings of P_m according to how c is covered,

D_s = D_L·D_R + ε·x·D_{L′}·D_R + (−ε)·x·D_L·D_{L′} = D_σ². ∎

## Lemma 2 (parity)

Over 𝔽₂ = ℤ/2, D_s ≡ f_m, and f_k(y + y²) = (1 + y)^{k+1} + y^{k+1} in 𝔽₂[y]. Consequently f_m ≡ 1 (mod 2) if and only if m + 1 is a power of 2.

Proof. Every s_e is odd, so D_s ≡ f_m. Over 𝔽₂, the roots of r² + r + (y + y²) are r = y and r = 1 + y, so g_k = (1 + y)^{k+1} + y^{k+1} satisfies g_k = g_{k−1} + (y + y²)·g_{k−2}, with g₀ = g₁ = 1; f_k satisfies the same recurrence with the same start. The substitution x ↦ y + y² is an injective ring homomorphism 𝔽₂[x] → 𝔽₂[y]. If m + 1 = 2^K, then (1 + y)^{2^K} = 1 + y^{2^K} and g_m = 1. If m + 1 is not a power of 2, some binomial coefficient C(m + 1, i) with 0 < i < m + 1 is odd, and g_m ≠ 1. ∎

This is Corless's dimension theorem.

## Lemma 3 (independence)

Let m = 2N + 1 be odd and c_e = x·f_{e−1}·f_{m−1−e} ∈ 𝔽₂[x]. Then c₁, …, c_N are linearly independent over 𝔽₂.

Proof. Drop the common factor x and substitute x = y + y². By Lemma 2,

f_{e−1}·f_{m−1−e} ↦ ((1 + y)^e + y^e)·((1 + y)^{m−e} + y^{m−e}) = A + φ_e,

with A = (1 + y)^m + y^m, whose constant term is 1, and φ_e = y^e(1 + y)^{m−e} + y^{m−e}(1 + y)^e. For 1 ≤ e ≤ N we have m − e ≥ N + 1 > e, so the lowest term of φ_e is y^e and φ_e has no constant term. Suppose Σ_{e∈S} c_e = 0 for a nonempty S ⊆ {1, …, N}. The constant term shows that |S| is even, so Σ_{e∈S} φ_e = 0. But the coefficient of y^{min S} in this sum is 1: it comes from φ_{min S} only, because every other φ_e has all its terms in degrees above min S. Contradiction. ∎

## Lemma 4 (mod 4 forces folding)

Let m = 2N + 1 = 2^K − 1 with K ≥ 2. If D_s ≡ 1 (mod 4), then s_N = −s_{N+1} and s_{m−e} = s_e for 1 ≤ e ≤ N − 1; that is, s = (σ, ε, −ε, rev σ).

Proof. Write s_e = 1 − 2t_e with t_e ∈ {0, 1}. Since ∏_{e∈M}(1 − 2t_e) ≡ 1 − 2·Σ_{e∈M} t_e (mod 4),

D_s ≡ f_m − 2·Σ_e t_e·x·f_{e−1}·f_{m−1−e}   (mod 4),

because the matchings of P_m containing edge e correspond to the matchings of the two paths left after deleting its endpoints, which have e − 1 and m − 1 − e vertices. By Lemma 2 all coefficients of f_m except the constant are even; write f_m = 1 + 2g. Then D_s ≡ 1 (mod 4) is equivalent to

Σ_{e=1}^{m−1} t_e·c_e = ḡ   in 𝔽₂[x],

where ḡ is g mod 2 and c_e is as in Lemma 3. Since c_e = c_{m−e}, the left side equals Σ_{e=1}^{N} w_e·c_e with w_e = t_e + t_{m−e} (mod 2). By Lemma 3 the vector (w₁, …, w_N) is uniquely determined by ḡ. Now every sequence of the form (σ, ε, −ε, rev σ) satisfies D_s ≡ 1 (mod 4): by Lemma 1, D_s = D_σ², and by Lemma 2, D_σ ≡ f_N ≡ 1 (mod 2) because N + 1 = 2^{K−1}, so D_σ = 1 + 2h and D_σ² = 1 + 4h + 4h² ≡ 1. Such a sequence has w_e = 0 for e < N and w_N = 1. Hence every s with D_s ≡ 1 (mod 4) has the same w, which is the claim. ∎

## Proof of the theorem

If T_m(u) is nilpotent, then D_s = 1, so D_s ≡ 1 (mod 2) and m = 2^K − 1 by Lemma 2. We show s ∈ F_K by induction on K. For K = 1, s is empty. For K ≥ 2, D_s = 1 implies D_s ≡ 1 (mod 4), so s = (σ, ε, −ε, rev σ) by Lemma 4, with σ of length N − 1 = 2^{K−1} − 2. By Lemma 1, D_σ² = D_s = 1, so (D_σ − 1)(D_σ + 1) = 0 in the integral domain ℤ[x]; as D_σ(0) = 1, D_σ = 1. Thus T_N(σ) is nilpotent, σ ∈ F_{K−1} by induction, and s ∈ F_K. Conversely, if s ∈ F_K then D_s = D_σ² = 1 by Lemma 1 and induction. ∎

## Remarks

- The exhaustive search at m = 63 is no longer needed once Lemma 4 is proved, but it extends Corless's computational verification from m = 31 to m = 63 independently of the proof.
- The argument never needs more than the reduction modulo 4. Modulo 4 the condition is already linear over 𝔽₂, and it alone forces the folded shape; the integer condition then passes to the half-size sequence through D_s = D_σ².
- D_s = 1 says that the Stieltjes continued fraction 1/(1 + s₁x/(1 + s₂x/(⋯/(1 + s_{2N}x)))) is a polynomial. Lemma 1 is analogous to the folding lemma for continued fractions (van der Poorten and Shallit, Folded continued fractions, J. Number Theory 40 (1992) 237–250); the exact correspondence has not been worked out here. For every member of F_K with K ≤ 6 that polynomial is 1 ± x ± x³ ± x⁷ ± ⋯ ± x^{2^{K−1}−1}, with all sign choices occurring (observed by computation, `sfraction.py`; not needed for the proof).
- Lemma 3 holds for every odd m. The power-of-two condition enters only through Lemma 2.
- The closed form in Lemma 2 is Binet's formula modulo 2. It is also known in the language of reversed Dickson polynomials of the second kind. It follows from their functional equation (Hong, Qin & Zhao, Finite Fields Appl. 37 (2016) 54–71), and Fang, Fernando & Wu (Involve 18 (2025) 297–328) state it for characteristic 2. Both were read in their arXiv versions. Lemmas 3 and 4 are what is new.

## Numerical checks

All scripts are in `computations/`. `computations/run_all.sh` reruns every check below in about half a minute on 14 cores; `computations/README.md` describes each script.

| Check | Script | Result |
|---|---|---|
| Brute force over all sign sequences, m ≤ 9 and m = 15 | `brute.py` | nilpotent sets equal F_K (2, 4, 8 at m = 3, 7, 15), none at other m |
| Exhaustive search of the mod-4 solution space, m = 7, 15, 31, 63 (2^31 candidates at m = 63; 14 cores, 20 s) | `gen_params.py`, `enum.c` | exactly 4, 8, 16, 32 nilpotent sequences, equal to F_K. Corless had checked up to m = 31 |
| Every mod-4 solution has the folded shape, m = 15 and 31 | `halves.py` | all 128 and 32,768 solutions are folded |
| Rank of the mod-4 system equals N, m ≤ 255 | `mod4.py` | full rank |
| Lemma 1 for 3,000 random σ of length ≤ 40 | `check_proof.py` | holds |
| Lemma 2 closed form, k < 200; parity of f_k, k < 4,096 | `check_proof.py` | holds |
| Lemma 3, m = 2^K − 1 up to 2,047 | `check_proof.py` | rank N |
| The matrices themselves, exact Gaussian-integer arithmetic: T^m = 0 and T^{m−1} ≠ 0 for every family member at m = 3, 7, 15; no other nilpotent at m = 7 (all 60) or among 200 random vectors at m = 15; det(λI − T) equals the recurrence for 60 random u | `check_matrices.py` | holds |

## Formal verification

The whole argument is formalised in Lean 4 (v4.36.0) with Mathlib (v4.36.0) in `lean/` (five files, 1,150 lines). This includes Corless's dimension theorem, the link from the matrices T_m(u) to the matching polynomial, and the count. The conjecture is stated in Corless's own terms. `CorlessVec k u` is his recursion, written literally on lists: the empty list, and u ↦ u ++ [1, i] ++ rev u or u ++ [i, 1] ++ rev u. His K is k + 1.

```lean
theorem Corless.corless_conjecture {u : List ℂ} (hu : ∀ x ∈ u, x = 1 ∨ x = Complex.I) :
    IsNilpotent (T (u.length + 1) (fun i => u.getD i 1)) ↔ ∃ k, CorlessVec k u

theorem Corless.card_nilpotent (k : ℕ) :
    {u : List ℂ | u.length = 2 ^ (k + 1) - 2 ∧ (∀ x ∈ u, x = 1 ∨ x = Complex.I) ∧
      IsNilpotent (T (u.length + 1) (fun i => u.getD i 1))}.ncard = 2 ^ k
```

`T m u` is the m × m matrix with zero diagonal, u at (i, i+1) and −u at (i+1, i). This is Corless's eq. (1), indexed from 0. The first theorem is Conjecture 5.1 together with Theorem 4.2. Since `CorlessVec k u` forces u to have length 2^(k+1) − 2, it also contains the dimension theorem. The second is the count 2^{K−1} in dimension 2^K − 1.

`lake build` (from a clean build directory) compiles all files with no errors and no warnings. Mathlib's style linters are on, except its copyright-header convention. `#print axioms` reports only Lean's standard axioms (propext, Classical.choice, Quot.sound) for both theorems, so there is no `sorry`.

The formal proof follows the lemmas above. The link to the matrices goes through Laplace expansion of the characteristic matrix: first along the first column, then along the first row of the second minor. This gives the identity charpoly(T_m(u)) = reflect_m(expand₂(D_s)). Mathlib's criterion then finishes the link: an endomorphism is nilpotent if and only if its characteristic polynomial is X^m. See `lean/README.md`.
