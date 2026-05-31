//! Hodge star operator ⋆ on differential forms.

use crate::{DifferentialForm, binomial, k_subsets};
use nalgebra::DMatrix;

/// Hodge star ⋆: Ω^k → Ω^(n-k).
#[derive(Clone, Debug)]
pub struct HodgeStar {
    pub manifold_dim: usize,
    pub from_degree: usize,
    matrix: DMatrix<f64>,
}

impl HodgeStar {
    pub fn new(manifold_dim: usize, from_degree: usize) -> Self {
        let matrix = compute_hodge_star_matrix(manifold_dim, from_degree);
        Self { manifold_dim, from_degree, matrix }
    }

    pub fn apply(&self, form: &DifferentialForm) -> DifferentialForm {
        let target_degree = self.manifold_dim - self.from_degree;
        let v = nalgebra::DVector::from_vec(form.coefficients.clone());
        let result = &self.matrix * v;
        DifferentialForm::new(self.manifold_dim, target_degree, result.iter().cloned().collect())
    }

    pub fn matrix(&self) -> &DMatrix<f64> { &self.matrix }

    pub fn verify_double_star(&self) -> bool {
        let k = self.from_degree;
        let n = self.manifold_dim;
        let star2 = Self::new(n, n - k);
        let double = star2.matrix() * &self.matrix;
        let expected_sign = if (k * (n - k)) % 2 == 0 { 1.0 } else { -1.0 };
        let dim = binomial(n, k);
        if dim == 0 { return true; }
        let identity = DMatrix::identity(dim, dim) * expected_sign;
        double.relative_eq(&identity, 1e-8, 1e-8)
    }
}

fn compute_hodge_star_matrix(n: usize, k: usize) -> DMatrix<f64> {
    let target_k = match n.checked_sub(k) {
        Some(t) => t,
        None => return DMatrix::zeros(0, binomial(n, k)),
    };
    let rows = binomial(n, target_k);
    let cols = binomial(n, k);
    if rows == 0 || cols == 0 {
        return DMatrix::zeros(rows, cols);
    }

    let basis_k = k_subsets(n, k);
    let basis_nk = k_subsets(n, target_k);
    let mut data = vec![0.0; rows * cols];

    for (j, sigma) in basis_k.iter().enumerate() {
        let complement: Vec<usize> = (0..n).filter(|x| !sigma.contains(x)).collect();
        if let Some(i) = basis_nk.iter().position(|t| *t == complement) {
            let mut perm: Vec<usize> = sigma.iter().cloned().collect();
            perm.extend(&complement);
            let sign = permutation_sign(&perm);
            data[i * cols + j] = sign as f64;
        }
    }
    DMatrix::from_vec(rows, cols, data)
}

fn permutation_sign(perm: &[usize]) -> i32 {
    let n = perm.len();
    let mut visited = vec![false; n];
    let mut sign = 1i32;
    for i in 0..n {
        if visited[i] { continue; }
        let mut cycle_len = 0;
        let mut j = i;
        while !visited[j] {
            visited[j] = true;
            j = perm[j];
            cycle_len += 1;
        }
        if cycle_len % 2 == 0 { sign *= -1; }
    }
    sign
}
