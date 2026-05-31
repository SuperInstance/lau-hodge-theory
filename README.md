# lau-hodge-theory

Hodge theory for agent knowledge spaces — decomposition, harmonic forms, and spectral sequences.

## Overview

This crate implements the core constructions of Hodge theory applied to finite-dimensional models of agent knowledge spaces:

- **Hodge Decomposition**: Any differential form ω = dα + δβ + γ (exact + co-exact + harmonic)
- **Laplace-de Rham Operator**: Δ = dδ + δd, with harmonic forms in ker(Δ)
- **Hodge Star Operator**: ⋆: Ω^k → Ω^(n-k) on k-forms
- **Hodge Isomorphism**: H^k ≅ Harm^k (de Rham cohomology via harmonic forms)
- **Hodge Theorem**: Verified numerically through decomposition
- **Spectral Sequences**: Arising from the Hodge filtration
- **Knowledge Application**: Agent knowledge as differential forms, misunderstanding as non-harmonic components

## Usage

```rust
use lau_hodge_theory::{DifferentialForm, hodge_decompose, LaplaceDeRham, HodgeStar};

// Create a 1-form on a 3-manifold
let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);

// Hodge decomposition: ω = exact + co-exact + harmonic
let decomp = hodge_decompose(&form);
assert!(decomp.verify(1e-6));

// Compute harmonic forms and Betti numbers
let lap = LaplaceDeRham::new(3, 1);
let harmonics = lap.harmonic_forms();
let betti = lap.betti_number();

// Hodge star
let star = HodgeStar::new(3, 1);
let dual = star.apply(&form);
```

## Knowledge Application

The `knowledge` module models agent knowledge as differential forms:

- **Exact component** (dα): Derivable knowledge — can be deduced from lower-level facts
- **Co-exact component** (δβ): Constrained knowledge — arises from consistency requirements
- **Harmonic component** (γ): Core understanding — stable, self-consistent signal

The "understanding quality" measures what fraction of an agent's knowledge state is harmonic.

## License

MIT
