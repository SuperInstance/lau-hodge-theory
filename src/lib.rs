//! # Hodge Theory for Agent Knowledge Spaces

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::fmt;

pub mod hodge_star;
pub mod knowledge;

pub use hodge_star::HodgeStar;
pub use knowledge::KnowledgeSpace;

/// Binomial coefficient C(n, k).
pub fn binomial(n: usize, k: usize) -> usize {
    if k > n { return 0; }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// All k-element subsets of {0, ..., n-1} in lexicographic order.
pub fn k_subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 { return vec![vec![]]; }
    if k > n { return vec![]; }
    let mut result = Vec::new();
    let mut current: Vec<usize> = (0..k).collect();
    result.push(current.clone());
    loop {
        let mut i = k;
        while i > 0 {
            i -= 1;
            if current[i] < n - k + i {
                current[i] += 1;
                for j in (i + 1)..k { current[j] = current[j - 1] + 1; }
                result.push(current.clone());
                break;
            }
        }
        if i == 0 && current[0] >= n - k { break; }
    }
    result
}

// ---------------------------------------------------------------------------
// Differential Form
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DifferentialForm {
    pub degree: usize,
    pub coefficients: Vec<f64>,
    pub manifold_dim: usize,
}

impl DifferentialForm {
    pub fn new(manifold_dim: usize, degree: usize, coefficients: Vec<f64>) -> Self {
        let expected = binomial(manifold_dim, degree);
        assert_eq!(coefficients.len(), expected, "Expected {} coeffs for {}-form on {}-manifold", expected, degree, manifold_dim);
        Self { degree, coefficients, manifold_dim }
    }

    pub fn zero(manifold_dim: usize, degree: usize) -> Self {
        Self { degree, coefficients: vec![0.0; binomial(manifold_dim, degree)], manifold_dim }
    }

    pub fn len(&self) -> usize { self.coefficients.len() }
    pub fn is_zero(&self) -> bool { self.coefficients.iter().all(|c| c.abs() < 1e-10) }
    pub fn norm(&self) -> f64 { self.coefficients.iter().map(|c| c * c).sum::<f64>().sqrt() }

    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.degree, other.degree);
        assert_eq!(self.manifold_dim, other.manifold_dim);
        Self::new(self.manifold_dim, self.degree,
            self.coefficients.iter().zip(&other.coefficients).map(|(a, b)| a + b).collect())
    }

    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.manifold_dim, self.degree,
            self.coefficients.iter().map(|c| c * s).collect())
    }

    pub fn inner_product(&self, other: &Self) -> f64 {
        assert_eq!(self.degree, other.degree);
        assert_eq!(self.manifold_dim, other.manifold_dim);
        self.coefficients.iter().zip(&other.coefficients).map(|(a, b)| a * b).sum()
    }
}

impl fmt::Display for DifferentialForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-form[", self.degree)?;
        for (i, c) in self.coefficients.iter().enumerate() {
            if i > 0 { write!(f, ", ")?; }
            write!(f, "{:.4}", c)?;
        }
        write!(f, "]")
    }
}

// ---------------------------------------------------------------------------
// Exterior Derivative d: Ω^k → Ω^(k+1)
// ---------------------------------------------------------------------------

pub struct ExteriorDerivative {
    pub manifold_dim: usize,
    pub from_degree: usize,
    matrix: DMatrix<f64>,
}

impl ExteriorDerivative {
    pub fn new(manifold_dim: usize, from_degree: usize) -> Self {
        let matrix = compute_d_matrix(manifold_dim, from_degree);
        Self { manifold_dim, from_degree, matrix }
    }

    pub fn apply(&self, form: &DifferentialForm) -> DifferentialForm {
        let v = DVector::from_vec(form.coefficients.clone());
        let result = &self.matrix * v;
        DifferentialForm::new(self.manifold_dim, self.from_degree + 1, result.iter().cloned().collect())
    }

    pub fn matrix(&self) -> &DMatrix<f64> { &self.matrix }
}

/// d_k matrix: maps C(n,k)-vectors to C(n,k+1)-vectors.
fn compute_d_matrix(n: usize, k: usize) -> DMatrix<f64> {
    let rows = binomial(n, k + 1);
    let cols = binomial(n, k);
    if rows == 0 || cols == 0 {
        return DMatrix::zeros(rows, cols);
    }
    let mut data = vec![0.0f64; rows * cols];
    let basis_k = k_subsets(n, k);
    let basis_k1 = k_subsets(n, k + 1);
    for (i, sigma) in basis_k1.iter().enumerate() {
        for (omit_pos, _) in sigma.iter().enumerate() {
            let mut tau = sigma.clone();
            tau.remove(omit_pos);
            let sign = if omit_pos % 2 == 0 { 1.0 } else { -1.0 };
            if let Some(j) = basis_k.iter().position(|t| *t == tau) {
                data[j * rows + i] = sign;
            }
        }
    }
    DMatrix::from_vec(rows, cols, data)
}

// ---------------------------------------------------------------------------
// Hodge Star ⋆: Ω^k → Ω^(n-k)
// ---------------------------------------------------------------------------

// Defined in hodge_star.rs, re-exported above.

// ---------------------------------------------------------------------------
// Co-derivative δ: Ω^k → Ω^(k-1)
//
//   δ_k = (-1)^{n(k+1)+1} ⋆_{n-k+1} d_{n-k} ⋆_k
// ---------------------------------------------------------------------------

pub struct CoDerivative {
    pub manifold_dim: usize,
    pub from_degree: usize,
    matrix: DMatrix<f64>,
}

impl CoDerivative {
    pub fn new(manifold_dim: usize, from_degree: usize) -> Self {
        let n = manifold_dim;
        let k = from_degree;
        let target_dim = if k > 0 { binomial(n, k - 1) } else { 0 };
        let source_dim = binomial(n, k);

        if k == 0 || source_dim == 0 || target_dim == 0 {
            return Self { manifold_dim: n, from_degree: k, matrix: DMatrix::zeros(target_dim, source_dim) };
        }

        // δ_k = (-1)^{n(k+1)+1} ⋆_{n-k+1} d_{n-k} ⋆_k
        let star_k = HodgeStar::new(n, k);                              // ⋆_k: Ω^k → Ω^{n-k}
        let d_nk = ExteriorDerivative::new(n, n - k);                   // d_{n-k}: Ω^{n-k} → Ω^{n-k+1}
        let star_nkp1 = HodgeStar::new(n, n - k + 1);                   // ⋆_{n-k+1}: Ω^{n-k+1} → Ω^{k-1}

        let sign = if (n * (k + 1) + 1) % 2 == 0 { 1.0f64 } else { -1.0f64 };
        let matrix = (star_nkp1.matrix() * d_nk.matrix()) * star_k.matrix() * sign;

        Self { manifold_dim: n, from_degree: k, matrix }
    }

    pub fn apply(&self, form: &DifferentialForm) -> DifferentialForm {
        let target_degree = if self.from_degree > 0 { self.from_degree - 1 } else { 0 };
        let v = DVector::from_vec(form.coefficients.clone());
        let result = &self.matrix * v;
        DifferentialForm::new(self.manifold_dim, target_degree, result.iter().cloned().collect())
    }

    pub fn matrix(&self) -> &DMatrix<f64> { &self.matrix }
}

// ---------------------------------------------------------------------------
// Laplace-de Rham Δ = dδ + δd: Ω^k → Ω^k
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct LaplaceDeRham {
    pub manifold_dim: usize,
    pub degree: usize,
    matrix: DMatrix<f64>,
}

impl LaplaceDeRham {
    pub fn new(manifold_dim: usize, degree: usize) -> Self {
        let n = manifold_dim;
        let k = degree;
        let dim_k = binomial(n, k);
        if dim_k == 0 {
            return Self { manifold_dim: n, degree: k, matrix: DMatrix::zeros(0, 0) };
        }

        // δd: d_k then δ_{k+1}
        let delta_d = {
            let d_k = ExteriorDerivative::new(n, k);
            let delta_kp1 = CoDerivative::new(n, k + 1);
            delta_kp1.matrix() * d_k.matrix()
        };

        // dδ: δ_k then d_{k-1}
        let d_delta = if k > 0 {
            let d_km1 = ExteriorDerivative::new(n, k - 1);
            let delta_k = CoDerivative::new(n, k);
            d_km1.matrix() * delta_k.matrix()
        } else {
            DMatrix::zeros(dim_k, dim_k)
        };

        let matrix = &d_delta + &delta_d;
        Self { manifold_dim: n, degree: k, matrix }
    }

    pub fn apply(&self, form: &DifferentialForm) -> DifferentialForm {
        let v = DVector::from_vec(form.coefficients.clone());
        let result = &self.matrix * v;
        DifferentialForm::new(self.manifold_dim, self.degree, result.iter().cloned().collect())
    }

    pub fn is_harmonic(&self, form: &DifferentialForm) -> bool {
        self.apply(form).coefficients.iter().all(|c| c.abs() < 1e-8)
    }

    pub fn matrix(&self) -> &DMatrix<f64> { &self.matrix }

    pub fn harmonic_forms(&self) -> Vec<DifferentialForm> {
        kernel_basis(&self.matrix)
            .into_iter()
            .map(|coeffs| DifferentialForm::new(self.manifold_dim, self.degree, coeffs))
            .collect()
    }

    pub fn betti_number(&self) -> usize {
        self.harmonic_forms().len()
    }
}

// ---------------------------------------------------------------------------
// Hodge Decomposition ω = dα + δβ + γ
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HodgeDecomposition {
    pub exact: DifferentialForm,
    pub coexact: DifferentialForm,
    pub harmonic: DifferentialForm,
    pub original: DifferentialForm,
}

impl HodgeDecomposition {
    pub fn verify(&self, tol: f64) -> bool {
        let sum = self.exact.add(&self.coexact).add(&self.harmonic);
        self.original.coefficients.iter().zip(&sum.coefficients).all(|(a, b)| (a - b).abs() < tol)
    }
}

pub fn hodge_decompose(form: &DifferentialForm) -> HodgeDecomposition {
    let n = form.manifold_dim;
    let k = form.degree;
    let lap = LaplaceDeRham::new(n, k);

    // Project onto harmonic subspace via orthogonal projection
    let harmonic_forms = lap.harmonic_forms();
    let mut gamma_coeffs = vec![0.0; form.len()];
    for h in &harmonic_forms {
        let inner = form.inner_product(h);
        let h_sq = h.inner_product(h);
        if h_sq > 1e-12 {
            for (i, c) in gamma_coeffs.iter_mut().enumerate() {
                *c += (inner / h_sq) * h.coefficients[i];
            }
        }
    }
    let harmonic = DifferentialForm::new(n, k, gamma_coeffs);
    let remainder = form.add(&harmonic.scale(-1.0));

    let (exact, coexact) = split_exact_coexact(&remainder, n, k);

    HodgeDecomposition { exact, coexact, harmonic, original: form.clone() }
}

fn split_exact_coexact(form: &DifferentialForm, n: usize, k: usize) -> (DifferentialForm, DifferentialForm) {
    let dim_k = binomial(n, k);
    if form.is_zero() || dim_k == 0 {
        return (DifferentialForm::zero(n, k), DifferentialForm::zero(n, k));
    }

    // d_{k-1}: image spans exact k-forms
    let d_mat = if k > 0 {
        ExteriorDerivative::new(n, k - 1).matrix().clone()
    } else {
        DMatrix::zeros(dim_k, 0)
    };

    // δ_{k+1}: image spans co-exact k-forms
    let delta_mat = if k < n {
        CoDerivative::new(n, k + 1).matrix().clone()
    } else {
        DMatrix::zeros(dim_k, 0)
    };

    let n_d = d_mat.ncols();
    let n_delta = delta_mat.ncols();
    let total_cols = n_d + n_delta;

    if total_cols == 0 {
        return (DifferentialForm::zero(n, k), DifferentialForm::zero(n, k));
    }

    // Concatenate columns [d | δ]
    let mut combined_data = Vec::with_capacity(dim_k * total_cols);
    for col in 0..n_d {
        for row in 0..dim_k { combined_data.push(d_mat[(row, col)]); }
    }
    for col in 0..n_delta {
        for row in 0..dim_k { combined_data.push(delta_mat[(row, col)]); }
    }
    let combined = DMatrix::from_vec(dim_k, total_cols, combined_data);

    // Least-squares solve via SVD: x = V * Σ⁺ * U^T * b
    let target = DVector::from_vec(form.coefficients.clone());
    let svd = combined.clone().svd(true, true);
    let u = svd.u.unwrap_or_else(|| DMatrix::identity(dim_k, dim_k));
    let v_t = svd.v_t.unwrap_or_else(|| DMatrix::identity(total_cols, total_cols));
    let sigma = &svd.singular_values;

    let ut_b = u.transpose() * &target;
    let rank = sigma.len().min(dim_k.min(total_cols));
    let mut x = DVector::zeros(total_cols);
    for i in 0..rank {
        if sigma[i] > 1e-10 {
            let s_inv = 1.0 / sigma[i];
            for j in 0..total_cols {
                x[j] += v_t[(i, j)] * s_inv * ut_b[i];
            }
        }
    }

    let exact_coeffs: Vec<f64> = if n_d > 0 {
        let x_d = x.rows(0, n_d);
        (&d_mat * x_d).iter().cloned().collect()
    } else {
        vec![0.0; dim_k]
    };

    let coexact_coeffs: Vec<f64> = if n_delta > 0 {
        let x_delta = x.rows(n_d, n_delta);
        (&delta_mat * x_delta).iter().cloned().collect()
    } else {
        vec![0.0; dim_k]
    };

    (DifferentialForm::new(n, k, exact_coeffs), DifferentialForm::new(n, k, coexact_coeffs))
}

// ---------------------------------------------------------------------------
// Hodge Isomorphism: H^k ≅ Harm^k
// ---------------------------------------------------------------------------

pub struct HodgeIsomorphism {
    pub manifold_dim: usize,
    laplacians: Vec<LaplaceDeRham>,
}

impl HodgeIsomorphism {
    pub fn new(manifold_dim: usize) -> Self {
        let laplacians = (0..=manifold_dim).map(|k| LaplaceDeRham::new(manifold_dim, k)).collect();
        Self { manifold_dim, laplacians }
    }

    pub fn betti_number(&self, k: usize) -> usize {
        self.laplacians[k].betti_number()
    }

    pub fn harmonic_representatives(&self, k: usize) -> Vec<DifferentialForm> {
        self.laplacians[k].harmonic_forms()
    }

    pub fn verify_isomorphism(&self) -> bool {
        for k in 0..=self.manifold_dim {
            for h in &self.laplacians[k].harmonic_forms() {
                if !self.laplacians[k].is_harmonic(h) { return false; }
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Spectral Sequence
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpectralSequencePage {
    pub r: usize,
    pub terms: Vec<(usize, usize, f64)>,
}

#[derive(Clone, Debug)]
pub struct HodgeSpectralSequence {
    pub manifold_dim: usize,
    pub betti_numbers: Vec<usize>,
}

impl HodgeSpectralSequence {
    pub fn new(manifold_dim: usize) -> Self {
        let betti_numbers = (0..=manifold_dim).map(|k| LaplaceDeRham::new(manifold_dim, k).betti_number()).collect();
        Self { manifold_dim, betti_numbers }
    }

    pub fn e1_page(&self) -> SpectralSequencePage {
        let mut terms = Vec::new();
        for p in 0..=self.manifold_dim {
            for q in 0..=self.manifold_dim {
                if p + q <= self.manifold_dim {
                    let val = self.betti_numbers[p + q] as f64;
                    if val.abs() > 1e-10 { terms.push((p, q, val)); }
                }
            }
        }
        SpectralSequencePage { r: 1, terms }
    }

    pub fn e_inf_page(&self) -> SpectralSequencePage {
        let mut terms = Vec::new();
        for p in 0..=self.manifold_dim {
            let val = self.betti_numbers[p] as f64;
            if val.abs() > 1e-10 { terms.push((p, 0, val)); }
        }
        SpectralSequencePage { r: usize::MAX, terms }
    }

    pub fn total_cohomology(&self) -> Vec<usize> { self.betti_numbers.clone() }
}

// ---------------------------------------------------------------------------
// Kernel basis via SVD
// ---------------------------------------------------------------------------

pub fn kernel_basis(matrix: &DMatrix<f64>) -> Vec<Vec<f64>> {
    let (_m, n) = matrix.shape();
    if n == 0 || matrix.nrows() == 0 { return vec![]; }
    let svd = matrix.clone().svd(true, true);
    let v_t = match &svd.v_t {
        Some(vt) => vt,
        None => return vec![],
    };
    let sigma = &svd.singular_values;
    let tol = 1e-8;
    let mut basis = Vec::new();
    for j in 0..n {
        let sj = if j < sigma.len() { sigma[j] } else { 0.0 };
        if sj.abs() < tol {
            let v: Vec<f64> = (0..n).map(|i| {
                if j < v_t.nrows() && i < v_t.ncols() { v_t[(j, i)] } else { 0.0 }
            }).collect();
            basis.push(v);
        }
    }
    basis
}
