# Lean formalisation of Corless's conjecture

This is a machine-checked proof, in Lean 4 with Mathlib, of Conjecture 5.1 of R. M. Corless, *Skew-symmetric tridiagonal Bohemian matrices*, Maple Transactions 1(2) (2021), doi:10.5206/mt.v1i2.14360.

The conjecture states that the recursively constructed vectors give all nilpotent skew-symmetric tridiagonal matrices with superdiagonal entries in {1, i}. The formalisation also proves the number of such matrices. The informal proof is in `../proof.md`.

## Main statements

```lean
-- CorlessNilpotent/Vectors.lean
inductive CorlessVec : ℕ → List ℂ → Prop
  | base : CorlessVec 0 []
  | step₁ {k : ℕ} {u : List ℂ} :
      CorlessVec k u → CorlessVec (k + 1) (u ++ [1, Complex.I] ++ u.reverse)
  | step₂ {k : ℕ} {u : List ℂ} :
      CorlessVec k u → CorlessVec (k + 1) (u ++ [Complex.I, 1] ++ u.reverse)

theorem Corless.corless_conjecture {u : List ℂ} (hu : ∀ x ∈ u, x = 1 ∨ x = Complex.I) :
    IsNilpotent (T (u.length + 1) (fun i => u.getD i 1)) ↔ ∃ k, CorlessVec k u

theorem Corless.card_nilpotent (k : ℕ) :
    {u : List ℂ | u.length = 2 ^ (k + 1) - 2 ∧ (∀ x ∈ u, x = 1 ∨ x = Complex.I) ∧
      IsNilpotent (T (u.length + 1) (fun i => u.getD i 1))}.ncard = 2 ^ k
```

- **`T m u`.** `T m u : Matrix (Fin m) (Fin m) ℂ` has zero diagonal, `u i` at position `(i, i+1)` and `-u i` at `(i+1, i)`. This is Corless's matrix `T_m(u₁, …, u_{m-1})` of his eq. (1), indexed from 0, so the list `u` is his vector.
- **`CorlessVec k u`.** This is the rule of his Theorem 4.2, written literally. It starts from the empty vector (the 1 × 1 zero matrix) and passes from `u` to `(u, 1, i, rev u)` or `(u, i, 1, rev u)`. Corless's `k` is our `k + 1`: `u` has length `2^(k+1) - 2`, and `T` has dimension `2^(k+1) - 1`.
- **`corless_conjecture`.** This is Conjecture 5.1 together with Theorem 4.2. Since every `CorlessVec k u` has length `2^(k+1) - 2`, it also contains Corless's dimension theorem (Theorem 4.1).
- **`card_nilpotent`.** This counts vectors. Since `u` is the superdiagonal of `T`, distinct vectors give distinct matrices, so it is the count of Corless's corollary: `2^(K-1)` nilpotent matrices in dimension `2^K - 1`.

The intermediate statements are:

```lean
-- CorlessNilpotent/Matrix.lean: the same, for any u : ℕ → ℂ with values in {1, i}
theorem Corless.T_nilpotent_iff {u : ℕ → ℂ} (hu : ∀ i, u i = 1 ∨ u i = Complex.I) {m : ℕ}
    (hm : 1 ≤ m) :
    IsNilpotent (T m u) ↔
      ∃ k, m = 2 ^ (k + 1) - 1 ∧ Fam k (fun i => if u i = 1 then 1 else -1)

-- CorlessNilpotent/Main.lean: the classification for the matching polynomial
theorem Corless.corless {s : ℕ → ℤ} (hs : IsSign s) {m : ℕ} (hm : 1 ≤ m) :
    Dv s m = 1 ↔ ∃ k, m = 2 ^ (k + 1) - 1 ∧ Fam k s
```

Here `Dv s m ∈ ℤ[X]` is the matching polynomial of the path with `m` vertices, whose edge `i` carries the sign `s i = u i²`. `Fam k s` is Corless's recursion on these signs. It requires two opposite signs in the middle, a right half that mirrors the left half, and a left half that is again in the family.

## Files

| File | Content |
|---|---|
| `CorlessNilpotent/Path.lean` | the matching polynomial `Dv`, its recurrences, the edge-split formula and reversal invariance |
| `CorlessNilpotent/Parity.lean` | Lemma 1 (folding identity `D = D_half²`), reduction mod 2, the closed form under `x = y + y²` over 𝔽₂, and Corless's dimension theorem |
| `CorlessNilpotent/Main.lean` | telescoping, the reduction mod 4, Lemma 3 (independence over 𝔽₂), Lemma 4 (mod 4 forces the folded shape) and the classification `corless` |
| `CorlessNilpotent/Matrix.lean` | the characteristic polynomial of `T m u` by Laplace expansion, its relation to `Dv`, and `T_nilpotent_iff` |
| `CorlessNilpotent/Vectors.lean` | Corless's vectors `CorlessVec`, their equivalence with `Fam`, `corless_conjecture` and `card_nilpotent` |
| `Axioms.lean` | prints the axioms used by the main theorems |

## Checking

You need [elan](https://github.com/leanprover/elan). The toolchain is `leanprover/lean4:v4.36.0` and Mathlib is `v4.36.0`.

```
lake exe cache get      # download the compiled Mathlib (about 5 GB)
lake build              # compiles the proof
lake env lean Axioms.lean
```

**Build.** `lake build` compiles all files with no errors and no warnings. The files use the module system, and Mathlib's standard linter set is on. The copyright-header linter is switched off in `lakefile.toml`, because it enforces a Mathlib convention.

**Axioms.** The last command prints `[propext, Classical.choice, Quot.sound]` for every main theorem. These are only Lean's standard axioms, so there is no `sorry`.

Checked on 17 November 2026.
