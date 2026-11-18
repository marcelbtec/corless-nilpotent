# Corless's nilpotent conjecture: proof, Lean formalisation and computations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24317787.svg)](https://doi.org/10.5281/zenodo.24317787) [![check](https://github.com/marcelbtec/corless-nilpotent/actions/workflows/check.yml/badge.svg)](https://github.com/marcelbtec/corless-nilpotent/actions/workflows/check.yml)

This repository proves Conjecture 5.1 of R. M. Corless, *Skew-symmetric tridiagonal Bohemian matrices*, Maple Transactions 1(2) (2021), Article 14360, [doi:10.5206/mt.v1i2.14360](https://doi.org/10.5206/mt.v1i2.14360). The proof is machine-checked in Lean 4 with Mathlib. The repository accompanies the paper *Nilpotent skew-symmetric tridiagonal Bohemian matrices: a proof of Corless's conjecture* by M. Blattner (arXiv preprint, 2026).

**Theorem.** Let T_m(u) be the m × m skew-symmetric tridiagonal matrix with zero diagonal, superdiagonal u = (u₁, …, u_{m−1}) and subdiagonal −u, where every u_k is 1 or i. Then T_m(u) is nilpotent if and only if two conditions hold:

- m = 2^K − 1 for some K ≥ 1;
- u is one of Corless's recursive vectors, built from the empty vector by u ↦ (u, 1, i, rev u) or u ↦ (u, i, 1, rev u).

There are exactly 2^(K−1) such matrices in dimension 2^K − 1, and none in other dimensions.

The key step reduces the problem modulo 4. There, a necessary condition for nilpotency is a linear system over 𝔽₂, and that system already forces the outermost step of the recursion. An exact identity, implicit in Corless's paper, then reduces the problem to half the dimension. `proof.md` has the full argument.

The paper contains further results: the classification of the nilpotent tridiagonal sign matrices (finite sections of the Feinberg–Zee hopping matrix), the (1,1) entry of the resolvent of a nilpotent T_m(u), and a second description of the nilpotent vectors by a doubling rule. These are not formalised. `computations/check_extras.py` checks them by exact computation in small cases.

## Contents

| Path | Content |
|---|---|
| `proof.md` | the informal proof (four lemmas and the theorem), the record of the openness check, and a table of numerical checks |
| `lean/` | the Lean 4 formalisation; `lean/README.md` lists the statements and files |
| `computations/` | Python and C scripts that check each lemma, the classification by brute force and by exhaustive search up to m = 63, the matrices themselves, and the further results of the paper in small cases; `computations/README.md` describes each script |
| `.github/workflows/check.yml` | continuous integration: builds the Lean proof, checks its axioms, and reruns the computations |

## Main formal statements

```lean
theorem Corless.corless_conjecture {u : List ℂ} (hu : ∀ x ∈ u, x = 1 ∨ x = Complex.I) :
    IsNilpotent (T (u.length + 1) (fun i => u.getD i 1)) ↔ ∃ k, CorlessVec k u

theorem Corless.card_nilpotent (k : ℕ) :
    {u : List ℂ | u.length = 2 ^ (k + 1) - 2 ∧ (∀ x ∈ u, x = 1 ∨ x = Complex.I) ∧
      IsNilpotent (T (u.length + 1) (fun i => u.getD i 1))}.ncard = 2 ^ k
```

`T` is Corless's matrix and `CorlessVec` is his recursion, written literally. Lean's k is Corless's K − 1. Both theorems depend only on Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`).

## Checking

**Lean.** You need [elan](https://github.com/leanprover/elan). The toolchain is `leanprover/lean4:v4.36.0` and Mathlib is `v4.36.0`.

```
cd lean
lake exe cache get          # downloads the compiled Mathlib, about 5 GB
lake build
lake env lean Axioms.lean   # prints the axioms of the main theorems
```

**Computations.** You need `python3`, which uses only the standard library, and a C compiler with 128-bit integers, such as gcc or clang. `make_figure.py` additionally needs matplotlib.

```
computations/run_all.sh     # about 30 s on 14 cores; set SLICES to the number of cores
```

## License

The code in `lean/` and `computations/` is licensed under the Apache License 2.0 (`LICENSE`), the license of Mathlib. The text of `proof.md`, the data files in `computations/` and the figure are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

See `CITATION.cff`. Each release is archived on Zenodo.

- All versions: [doi:10.5281/zenodo.24317787](https://doi.org/10.5281/zenodo.24317787)
- Version 1.0.0: [doi:10.5281/zenodo.24317788](https://doi.org/10.5281/zenodo.24317788)

> Blattner, M. (2026). *Corless's nilpotent conjecture: proof, Lean formalisation and computations* (Version 1.1.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.24317787
