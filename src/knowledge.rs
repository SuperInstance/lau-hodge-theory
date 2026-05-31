//! Application: Agent Knowledge as Differential Forms.
//!
//! Models agent knowledge states as differential forms on a knowledge manifold.
//! Misunderstanding corresponds to non-harmonic components.

use crate::{
    DifferentialForm, HodgeIsomorphism,
    hodge_decompose,
};
use serde::{Deserialize, Serialize};

/// A knowledge form: a differential form representing an agent's knowledge state.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KnowledgeForm {
    /// The underlying differential form.
    pub form: DifferentialForm,
    /// Semantic labels for the basis elements.
    pub labels: Vec<String>,
}

impl KnowledgeForm {
    /// Create a knowledge form from coefficients and labels.
    pub fn new(form: DifferentialForm, labels: Vec<String>) -> Self {
        assert_eq!(form.len(), labels.len());
        Self { form, labels }
    }

    /// Decompose knowledge into:
    /// - Exact: derivable knowledge (can be deduced from lower-level facts)
    /// - Co-exact: constrained knowledge (arises from consistency requirements)
    /// - Harmonic: core understanding (stable, self-consistent)
    pub fn decompose(&self) -> KnowledgeDecomposition {
        let hodge = hodge_decompose(&self.form);
        KnowledgeDecomposition {
            derivable: hodge.exact,
            constrained: hodge.coexact,
            core_understanding: hodge.harmonic,
            total: self.form.clone(),
        }
    }

    /// Compute the "understanding quality" — fraction of knowledge that is harmonic.
    pub fn understanding_quality(&self) -> f64 {
        let decomp = self.decompose();
        let total_norm = self.form.norm();
        if total_norm < 1e-10 {
            return 1.0;
        }
        decomp.core_understanding.norm() / total_norm
    }
}

/// Decomposed knowledge state.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KnowledgeDecomposition {
    /// Derivable knowledge (exact component dα).
    pub derivable: DifferentialForm,
    /// Constrained knowledge (co-exact component δβ).
    pub constrained: DifferentialForm,
    /// Core understanding (harmonic component γ).
    pub core_understanding: DifferentialForm,
    /// Total knowledge form.
    pub total: DifferentialForm,
}

impl KnowledgeDecomposition {
    /// Verify the decomposition.
    pub fn verify(&self) -> bool {
        let sum = self.derivable.add(&self.constrained).add(&self.core_understanding);
        self.total
            .coefficients
            .iter()
            .zip(&sum.coefficients)
            .all(|(a, b)| (a - b).abs() < 1e-6)
    }
}

/// A knowledge space: the full space of knowledge forms on a topic manifold.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KnowledgeSpace {
    /// Dimension of the knowledge manifold.
    pub manifold_dim: usize,
    /// Betti numbers.
    pub betti_numbers: Vec<usize>,
}

impl KnowledgeSpace {
    /// Create a knowledge space of given dimension.
    pub fn new(manifold_dim: usize) -> Self {
        let iso = HodgeIsomorphism::new(manifold_dim);
        let betti_numbers = (0..=manifold_dim).map(|k| iso.betti_number(k)).collect();
        Self {
            manifold_dim,
            betti_numbers,
        }
    }

    /// The Euler characteristic: χ = Σ (-1)^k b_k.
    pub fn euler_characteristic(&self) -> i32 {
        self.betti_numbers
            .iter()
            .enumerate()
            .map(|(k, b)| if k % 2 == 0 { *b as i32 } else { -(*b as i32) })
            .sum()
    }

    /// Create a zero knowledge form of given degree.
    pub fn zero_form(&self, degree: usize) -> KnowledgeForm {
        KnowledgeForm {
            form: DifferentialForm::zero(self.manifold_dim, degree),
            labels: vec![String::new(); crate::binomial(self.manifold_dim, degree)],
        }
    }

    /// Analyze a knowledge form: decompose and measure quality.
    pub fn analyze(&self, form: &KnowledgeForm) -> KnowledgeAnalysis {
        let decomp = form.decompose();
        let quality = form.understanding_quality();
        KnowledgeAnalysis {
            decomposition: decomp,
            understanding_quality: quality,
        }
    }
}

/// Analysis result for a knowledge form.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KnowledgeAnalysis {
    pub decomposition: KnowledgeDecomposition,
    pub understanding_quality: f64,
}
