# lau-hodge-theory

Any differential form on a compact manifold decomposes into three orthogonal pieces: exact, co-exact, and harmonic. That's the Hodge decomposition — ω = dα + δβ + γ — and the harmonic piece γ encodes the topology.

Harmonic forms are the notes that persist — the ones a hermit crab hums walking between tide pools.

## The math in 60 seconds

The **Laplace-de Rham operator** Δ = dδ + δd generalizes the Laplacian to differential forms. Its kernel — the **harmonic forms** — are the forms that are simultaneously closed (dω = 0) and co-closed (δω = 0).

Key results:

- **Hodge decomposition:** Ωᵏ = im(d) ⊕ im(δ) ⊕ ker(Δ) — every form splits uniquely
- **Hodge isomorphism:** Hᵏ(M) ≅ Harmᵏ(M) — cohomology equals harmonic forms
- **Hodge star:** ⋆: Ωᵏ → Ωⁿ⁻ᵏ, with ⋆² = (-1)ᵏ⁽ⁿ⁻ᵏ⁾
- **Spectral sequences:** pages E₁, E₂, ... converging to the cohomology from the Hodge filtration
- **Knowledge application:** agent understanding as a differential form — misunderstanding is the non-harmonic component

References: Warner, *Foundations of Differentiable Manifolds and Lie Groups* (1983)

## Quick start

```rust
use lau_hodge_theory::{HodgeDecomposition, LaplaceDeRham, HodgeStar};

// Create a Hodge decomposition on a 3-manifold
let hodge = HodgeDecomposition::new(3); // 3-dimensional manifold

// Decompose a 1-form into exact + co-exact + harmonic
let omega = vec![1.0, 2.0, 0.5, -1.0, 0.3, 1.5]; // 6-component 1-form on M³
let (exact, coexact, harmonic) = hodge.decompose(&omega);

// Verify orthogonality: exact · coexact = 0, etc.
assert!(hodge.verify_orthogonality(&exact, &coexact, &harmonic));

// Hodge star on 1-forms in 3D: ⋆: Ω¹ → Ω²
let star = HodgeStar::new(3);
let star_omega = star.apply(&omega, 1); // 1-form → 2-form

// Betti numbers from harmonic form count
let betti = hodge.betti_numbers(); // [b₀, b₁, b₂, b₃]
```

## Key types

| Type | What it is |
|------|-----------|
| `HodgeDecomposition` | Splits any form into exact + co-exact + harmonic via SVD |
| `LaplaceDeRham` | The operator Δ = dδ + δd with kernel computation |
| `HodgeStar` | The ⋆ operator mapping Ωᵏ → Ωⁿ⁻ᵏ |
| `HarmonicForms` | Basis for ker(Δ) — computes Betti numbers |
| `SpectralSequence` | E₁ page from Hodge filtration, converging to H*(M) |

## Contributing

[Open an issue](https://github.com/SuperInstance/lau-hodge-theory/issues) or PR. Interesting directions:

- Non-compact manifolds (weighted L² spaces)
- Hodge theory on simplicial complexes (discrete exterior calculus)
- Applications to data analysis (the "Hodge filter")
