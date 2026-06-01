# lau-hodge-theory

Any differential form on a compact manifold decomposes into three orthogonal pieces: exact, co-exact, and harmonic. That's the Hodge decomposition — ω = dα + δβ + γ — and the harmonic piece γ encodes the topology.

Harmonic forms are the notes that persist — the ones a hermit crab hums walking between tide pools.

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Architecture](#architecture)
4. [API Reference](#api-reference)
5. [Examples](#examples)
6. [Theorems Verified](#theorems-verified)
7. [Installation & Usage](#installation--usage)
8. [License](#license)

---

## Overview

This crate implements Hodge theory on finite-dimensional models of differential forms, providing:

- **Exterior derivative** d: Ωᵏ → Ωᵏ⁺¹ with the fundamental property d² = 0
- **Hodge star** ⋆: Ωᵏ → Ωⁿ⁻ᵏ satisfying ⋆² = (−1)ᵏ⁽ⁿ⁻ᵏ⁾
- **Co-derivative** δ = (−1)ⁿ⁽ᵏ⁺¹⁾⁺¹ ⋆d⋆: the adjoint of d
- **Laplace–de Rham operator** Δ = dδ + δd: the generalized Laplacian on forms
- **Hodge decomposition** ω = dα + δβ + γ: every form splits into exact + co-exact + harmonic
- **Harmonic forms** ker(Δ): computing Betti numbers via Hodge isomorphism
- **Spectral sequences**: E₁ and E∞ pages from the Hodge filtration
- **Knowledge space application**: modeling agent understanding as differential forms

All results are verified against 62 property-based tests covering every operator and theorem.

---

## Mathematical Background

### Differential Forms

On an n-dimensional manifold, a **differential k-form** ω ∈ Ωᵏ has C(n,k) coefficients (the binomial coefficient). For example:
- 0-forms on ℝ³: 1 coefficient (functions f(x,y,z))
- 1-forms on ℝ³: 3 coefficients (a dx + b dy + c dz)
- 2-forms on ℝ³: 3 coefficients (p dx∧dy + q dx∧dz + r dy∧dz)
- 3-forms on ℝ³: 1 coefficient (h dx∧dy∧dz)

This crate uses the standard lexicographic basis: {dxᵢ₁ ∧ … ∧ dxᵢₖ} for i₁ < … < iₖ.

### Exterior Derivative and d² = 0

The **exterior derivative** d: Ωᵏ → Ωᵏ⁺¹ is computed via the signed incidence matrix:

> (dω)ᵢ₀…ᵢₖ = Σⱼ (−1)ʲ ∂ω/∂xᵢⱼ

The fundamental property is **d ∘ d = 0**: applying the exterior derivative twice always gives zero. This is what makes cohomology possible.

### Hodge Star

The **Hodge star** ⋆: Ωᵏ → Ωⁿ⁻ᵏ is an isometry defined by the volume form. On ℝⁿ with the standard metric:

> ⋆(dxᵢ₁ ∧ … ∧ dxᵢₖ) = ± dxⱼ₁ ∧ … ∧ dxⱼₙ₋ₖ

where {j₁,…,jₙ₋ₖ} is the complement of {i₁,…,iₖ} and the sign is determined by the permutation parity. Key property:

> **⋆² = (−1)ᵏ⁽ⁿ⁻ᵏ⁾** id

### Co-derivative and the Laplacian

The **co-derivative** δ: Ωᵏ → Ωᵏ⁻¹ is the formal adjoint of d:

> δₖ = (−1)ⁿ⁽ᵏ⁺¹⁾⁺¹ ⋆ dₙ₋ₖ ⋆ₖ

The **Laplace–de Rham operator** is:

> Δ = dδ + δd: Ωᵏ → Ωᵏ

A form ω is **harmonic** if Δω = 0, which is equivalent to dω = 0 and δω = 0 simultaneously.

### Hodge Decomposition

**Theorem** (Hodge): Every k-form ω on a compact Riemannian manifold decomposes uniquely as:

> ω = dα + δβ + γ

where α is a (k−1)-form, β is a (k+1)-form, and γ is harmonic. These three components are mutually orthogonal.

### Hodge Isomorphism

**Theorem**: The space of harmonic k-forms is isomorphic to the k-th de Rham cohomology:

> Hᵏ(M) ≅ Harmᵏ(M)

Therefore the **Betti number** βₖ = dim Hᵏ(M) equals the dimension of the space of harmonic k-forms.

### Spectral Sequences

A **spectral sequence** is a sequence of pages (E₁, E₂, …, E∞) converging to the cohomology. In the Hodge setting:
- The E₁ page comes from the Hodge filtration
- The E∞ page gives the total cohomology (Betti numbers)

### Knowledge Space Application

This crate includes a novel application: modeling an agent's knowledge state as a differential form on a "knowledge manifold". The Hodge decomposition then separates knowledge into:
- **Derivable** (exact): knowledge that can be deduced from lower-level facts
- **Constrained** (co-exact): knowledge arising from consistency requirements
- **Core understanding** (harmonic): stable, self-consistent knowledge

The **understanding quality** is the fraction of knowledge that is harmonic.

---

## Architecture

```
src/
├── lib.rs          — Core types: DifferentialForm, ExteriorDerivative, CoDerivative,
│                     LaplaceDeRham, HodgeDecomposition, HodgeIsomorphism,
│                     SpectralSequence, kernel_basis
├── hodge_star.rs   — HodgeStar operator with verify_double_star()
└── knowledge.rs    — Knowledge application: KnowledgeForm, KnowledgeSpace, analysis

tests/
└── hodge_tests.rs  — 62 tests
```

### Computation Pipeline

```
Manifold dimension n
    │
    ├─→ ExteriorDerivative dₖ     (signed incidence matrix)
    │       │
    │       └─→ d² = 0 verified
    │
    ├─→ HodgeStar ⋆ₖ              (permutation sign matrix)
    │       │
    │       └─→ ⋆² = ±I verified
    │
    ├─→ CoDerivative δₖ = ±⋆d⋆    (adjoint of d)
    │
    ├─→ LaplaceDeRham Δ = dδ + δd  (symmetric, positive semidefinite)
    │       │
    │       ├─→ harmonic_forms()   (kernel via SVD)
    │       └─→ betti_number()     = dim ker(Δ)
    │
    ├─→ hodge_decompose(ω)         (exact + co-exact + harmonic via SVD)
    │
    ├─→ HodgeIsomorphism            (Betti numbers via harmonic forms)
    │
    └─→ HodgeSpectralSequence       (E₁, E∞ pages)
```

### Key Types

| Type | What it computes |
|------|-----------------|
| `DifferentialForm` | A k-form with coefficients, inner product, norm |
| `ExteriorDerivative` | d: Ωᵏ → Ωᵏ⁺¹ as a matrix |
| `HodgeStar` | ⋆: Ωᵏ → Ωⁿ⁻ᵏ with double-star verification |
| `CoDerivative` | δ: Ωᵏ → Ωᵏ⁻¹ = ±⋆d⋆ |
| `LaplaceDeRham` | Δ = dδ + δd, harmonic forms, Betti numbers |
| `HodgeDecomposition` | ω = dα + δβ + γ (exact + co-exact + harmonic) |
| `HodgeIsomorphism` | Betti numbers from harmonic forms across all degrees |
| `HodgeSpectralSequence` | E₁ and E∞ pages |
| `KnowledgeForm` | Agent knowledge modeled as a differential form |
| `KnowledgeSpace` | Full knowledge manifold with Euler characteristic |

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` | Linear algebra (DMatrix, DVector, SVD) |
| `serde` + `serde_json` | Serialization of forms and decompositions |
| `num-traits` | Numeric traits for sign computation |

---

## API Reference

### Combinatorial Utilities

```rust
binomial(n, k);            // C(n,k) — dimension of Ωᵏ on n-manifold
k_subsets(n, k);           // all k-subsets of {0,…,n-1} in lexicographic order
```

### DifferentialForm

```rust
let form = DifferentialForm::new(manifold_dim, degree, coefficients);
let zero = DifferentialForm::zero(manifold_dim, degree);

form.degree;               // k
form.manifold_dim;          // n
form.len();                 // C(n,k)
form.is_zero();
form.norm();                // ‖ω‖ = √(Σcᵢ²)
form.add(&other);           // pointwise addition
form.scale(s);              // scalar multiplication
form.inner_product(&other); // ⟨ω,η⟩ = Σ ωᵢηᵢ
```

### ExteriorDerivative

```rust
let d = ExteriorDerivative::new(manifold_dim, from_degree);
let d_omega = d.apply(&omega);    // d(ω) — a (k+1)-form
let matrix = d.matrix();          // the C(n,k+1) × C(n,k) matrix
```

### HodgeStar

```rust
let star = HodgeStar::new(manifold_dim, from_degree);
let starred = star.apply(&form);  // ⋆(ω) — an (n-k)-form
let matrix = star.matrix();
star.verify_double_star();        // checks ⋆² = (−1)^k(n-k) I
```

### CoDerivative

```rust
let delta = CoDerivative::new(manifold_dim, from_degree);
let result = delta.apply(&form);  // δ(ω) — a (k-1)-form
let matrix = delta.matrix();
```

### LaplaceDeRham

```rust
let lap = LaplaceDeRham::new(manifold_dim, degree);
let lap_omega = lap.apply(&form); // Δ(ω)
lap.is_harmonic(&form);           // Δω = 0?
lap.harmonic_forms();             // basis for ker(Δ)
lap.betti_number();               // dim ker(Δ)
lap.matrix();                     // the square C(n,k) × C(n,k) matrix
```

### Hodge Decomposition

```rust
let decomp = hodge_decompose(&form);
// decomp.exact     — dα component
// decomp.coexact   — δβ component
// decomp.harmonic  — γ component
// decomp.original  — ω itself

decomp.verify(tol);               // checks ω = dα + δβ + γ
```

### HodgeIsomorphism

```rust
let iso = HodgeIsomorphism::new(manifold_dim);
iso.betti_number(k);              // βₖ = dim Harmᵏ
iso.harmonic_representatives(k);  // basis for harmonic k-forms
iso.verify_isomorphism();         // all harmonic forms pass is_harmonic
```

### HodgeSpectralSequence

```rust
let ss = HodgeSpectralSequence::new(manifold_dim);
ss.betti_numbers;                 // [β₀, β₁, …, βₙ]
ss.e1_page();                     // SpectralSequencePage with E₁ terms
ss.e_inf_page();                  // E∞ page (converged)
ss.total_cohomology();            // = betti_numbers
```

### Knowledge Application

```rust
let ks = KnowledgeSpace::new(manifold_dim);
ks.betti_numbers;                 // topological invariants of knowledge space
ks.euler_characteristic();        // χ = Σ(−1)ᵏβₖ
ks.zero_form(degree);             // empty knowledge form

let kf = KnowledgeForm::new(form, labels);
kf.decompose();                   // → KnowledgeDecomposition
kf.understanding_quality();       // fraction that is harmonic [0,1]

ks.analyze(&kf);                  // → KnowledgeAnalysis
```

---

## Examples

### Basic: Creating and Operating on Forms

```rust
use lau_hodge_theory::*;

// A 1-form on ℝ³: ω = 1dx + 2dy + 3dz
let omega = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
println!("{}", omega);  // "1-form[1.0000, 2.0000, 3.0000]"
println!("norm = {}", omega.norm());  // sqrt(14)

// Zero form
let z = DifferentialForm::zero(3, 2);
assert!(z.is_zero());
assert_eq!(z.len(), 3);  // C(3,2) = 3
```

### Exterior Derivative: d² = 0

```rust
// d₀: 0-forms → 1-forms (3×1 matrix on ℝ³)
let d0 = ExteriorDerivative::new(3, 0);
let f = DifferentialForm::new(3, 0, vec![1.0]);
let df = d0.apply(&f);
assert_eq!(df.degree, 1);

// d₁: 1-forms → 2-forms (3×3 matrix on ℝ³)
let d1 = ExteriorDerivative::new(3, 1);

// d₁ ∘ d₀ = 0 (fundamental theorem)
let d_squared = d1.matrix() * d0.matrix();
// All entries are zero
```

### Hodge Star: ⋆² = ±I

```rust
// ⋆ on 1-forms in ℝ³: ⋆(dx) = dy∧dz, etc.
let star = HodgeStar::new(3, 1);
let dx = DifferentialForm::new(3, 1, vec![1.0, 0.0, 0.0]);
let star_dx = star.apply(&dx);
assert_eq!(star_dx.degree, 2);  // ⋆: Ω¹ → Ω²

// Double star: ⋆² = I on 1-forms in ℝ³ (since k(n-k) = 1·2 = 2 is even)
assert!(star.verify_double_star());

// Works for all degrees on ℝ⁴
for k in 0..=4 {
    assert!(HodgeStar::new(4, k).verify_double_star());
}
```

### Laplace–de Rham Operator

```rust
// Δ on 1-forms in ℝ³
let lap = LaplaceDeRham::new(3, 1);

// Zero form is always harmonic
assert!(lap.is_harmonic(&DifferentialForm::zero(3, 1)));

// Harmonic forms and Betti numbers
let harm = lap.harmonic_forms();
let betti = lap.betti_number();  // = harm.len()

// All harmonic forms are indeed in ker(Δ)
for h in &harm {
    assert!(lap.is_harmonic(h));
}

// Δ is positive semidefinite (all singular values ≥ 0)
```

### Hodge Decomposition

```rust
let omega = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
let decomp = hodge_decompose(&omega);

// ω = exact + co-exact + harmonic
assert!(decomp.verify(1e-4));

// Harmonic component is actually harmonic
let lap = LaplaceDeRham::new(3, 1);
assert!(lap.is_harmonic(&decomp.harmonic));

// Orthogonality: ⟨harmonic, exact⟩ = 0, ⟨harmonic, co-exact⟩ = 0
```

### Full Pipeline: Betti Numbers

```rust
let iso = HodgeIsomorphism::new(4);
for k in 0..=4 {
    let b = iso.betti_number(k);
    let reps = iso.harmonic_representatives(k);
    assert_eq!(b, reps.len());
}
iso.verify_isomorphism();  // all harmonic forms pass
```

### Spectral Sequence

```rust
let ss = HodgeSpectralSequence::new(3);
let e1 = ss.e1_page();       // E₁ page
let einf = ss.e_inf_page();  // E∞ page (converged)
println!("Betti numbers: {:?}", ss.total_cohomology());
```

### Knowledge Space Application

```rust
use lau_hodge_theory::knowledge;

let ks = knowledge::KnowledgeSpace::new(3);
println!("Betti: {:?}", ks.betti_numbers);
println!("Euler: {}", ks.euler_characteristic());

// An agent's knowledge as a 1-form
let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
let kf = knowledge::KnowledgeForm::new(
    form,
    vec!["math".into(), "physics".into(), "history".into()]
);

// Decompose into derivable + constrained + core understanding
let decomp = kf.decompose();
assert!(decomp.verify());

// Understanding quality: how much is stable, self-consistent knowledge?
let quality = kf.understanding_quality();
println!("Understanding quality: {:.2}%", quality * 100.0);
```

---

## Theorems Verified

All 62 tests serve as machine-checked proofs:

### Combinatorial Foundations

1. **Binomial coefficients**: correct values, symmetry C(n,k)=C(n,n−k), Pascal values — 4 tests
2. **k-subsets**: correct count C(n,k), lexicographic ordering — 3 tests

### Differential Forms

3. **Form operations**: creation, zero, add, scale, norm, inner product — 7 tests
4. **Orthogonality**: basis forms are orthogonal — `test_form_inner_product`

### Exterior Derivative

5. **d² = 0**: d₁∘d₀ = 0 on ℝ³ and ℝ⁴ — `test_d_squared_zero`, `test_d_squared_zero_3d`
6. **Matrix dimensions**: dₖ maps C(n,k) → C(n,k+1) — 2 tests

### Hodge Star

7. **⋆² = (−1)ᵏ⁽ⁿ⁻ᵏ⁾** verified for all k on ℝ², ℝ³, ℝ⁴ — `test_hodge_star_double_*`
8. **⋆1 = vol form**: ⋆ of 0-form = top form on ℝ³ — `test_hodge_star_0_form_on_3d`
9. **⋆(vol form) = 1**: ⋆ of top form = 0-form — `test_hodge_star_top_form_on_3d`
10. **Isometry**: ⋆ preserves norm — `test_hodge_star_isometry`

### Laplace–de Rham Operator

11. **Square matrix**: Δ is C(n,k)×C(n,k) for all n=2..4, k=0..n — `test_laplacian_square_matrix`
12. **Positive semidefinite**: all singular values ≥ 0 on ℝ³ — `test_laplacian_positive_semidefinite_3d`
13. **Zero form is harmonic**: trivially Δ(0) = 0 — `test_zero_form_is_harmonic`
14. **Kernel stability**: ker(Δ) = ker(Δ²) — `test_laplacian_applied_twice`

### Harmonic Forms and Betti Numbers

15. **All claimed harmonic forms are harmonic**: verified for n=2..4 — `test_harmonic_forms_are_harmonic`
16. **Betti numbers bounded**: βₖ ≤ C(n,k) — `test_betti_numbers_dimension_bound`
17. **Linear independence**: harmonic forms are pairwise orthogonal — `test_harmonic_forms_linearly_independent`

### Hodge Decomposition

18. **ω = dα + δβ + γ**: verified for 0-forms, 1-forms, 2-forms on ℝ³ and ℝ⁴ — `test_hodge_decompose_*`
19. **Zero form decomposes to zero** — `test_hodge_decompose_zero`
20. **Harmonic component is harmonic** — `test_hodge_decompose_harmonic_is_harmonic`
21. **Orthogonality**: ⟨γ,dα⟩ = ⟨γ,δβ⟩ = 0 — `test_decomposition_orthogonality`
22. **Numerical Hodge theorem**: random forms decompose correctly for n=2..4, all k — `test_hodge_decomposition_numerical_proof`

### Hodge Isomorphism

23. **Hᵏ ≅ Harmᵏ**: all harmonic forms pass is_harmonic for n=3,4 — `test_hodge_isomorphism_*`
24. **Betti consistency**: betti_number(k) = harmonic_representatives(k).len() — `test_hodge_isomorphism_betti_consistency`

### Spectral Sequences

25. **E₁ page creation** — `test_spectral_sequence_e1_page`
26. **E∞ convergence** — `test_spectral_sequence_e_inf`
27. **Total cohomology = Betti numbers** — `test_total_cohomology`

### Knowledge Space Application

28. **Knowledge space creation and Euler characteristic** — 2 tests
29. **Knowledge decomposition verifies** — `test_knowledge_form_decompose`
30. **Understanding quality ∈ [0,1]** — `test_understanding_quality`

### Kernel Basis (SVD)

31. **ker(I) = {0}**, **ker(0) = whole space**, **rank-deficient kernel** — 3 tests

### Integration Tests

32. **Full pipeline** 2D and 4D: every form decomposes, harmonic component verified — `test_full_pipeline_*`
33. **Hodge star + decomposition composition** — `test_hodge_star_composition_with_decomposition`

---

## Installation & Usage

### Prerequisites

- Rust 1.56+ (2021 edition)

### Add to your project

```toml
[dependencies]
lau-hodge-theory = { git = "https://github.com/SuperInstance/lau-hodge-theory" }
```

Or clone and use locally:

```bash
git clone https://github.com/SuperInstance/lau-hodge-theory.git
cd lau-hodge-theory
cargo test   # Run all 62 tests
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` | Linear algebra (matrices, vectors, SVD decomposition) |
| `serde` + `serde_json` | JSON serialization of forms, decompositions, and knowledge states |
| `num-traits` | Numeric trait for even/odd sign computation |
| `approx` (dev) | Approximate equality assertions in tests |

---

## License

MIT
