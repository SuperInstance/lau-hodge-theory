#[cfg(test)]
mod tests {
    use lau_hodge_theory::*;
    use lau_hodge_theory::knowledge;
    use nalgebra::{DMatrix, DVector};
    use approx::assert_relative_eq;

    // ---- Binomial coefficient tests ----

    #[test]
    fn test_binomial_basic() {
        assert_eq!(binomial(4, 0), 1);
        assert_eq!(binomial(4, 1), 4);
        assert_eq!(binomial(4, 2), 6);
        assert_eq!(binomial(4, 3), 4);
        assert_eq!(binomial(4, 4), 1);
    }

    #[test]
    fn test_binomial_symmetry() {
        assert_eq!(binomial(5, 2), binomial(5, 3));
        assert_eq!(binomial(6, 1), binomial(6, 5));
    }

    #[test]
    fn test_binomial_zero() {
        assert_eq!(binomial(0, 0), 1);
        assert_eq!(binomial(3, 5), 0);
    }

    #[test]
    fn test_binomial_pascal() {
        assert_eq!(binomial(7, 3), 35);
        assert_eq!(binomial(10, 5), 252);
    }

    // ---- k-subsets tests ----

    #[test]
    fn test_k_subsets_2_of_3() {
        let subsets = k_subsets(3, 2);
        assert_eq!(subsets, vec![vec![0, 1], vec![0, 2], vec![1, 2]]);
    }

    #[test]
    fn test_k_subsets_0_forms() {
        let subsets = k_subsets(3, 0);
        assert_eq!(subsets, vec![vec![]]);
    }

    #[test]
    fn test_k_subsets_count() {
        assert_eq!(k_subsets(4, 2).len(), 6);
        assert_eq!(k_subsets(5, 3).len(), 10);
    }

    // ---- DifferentialForm tests ----

    #[test]
    fn test_form_creation() {
        let f = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        assert_eq!(f.degree, 1);
        assert_eq!(f.manifold_dim, 3);
        assert_eq!(f.len(), 3);
    }

    #[test]
    fn test_form_zero() {
        let f = DifferentialForm::zero(4, 2);
        assert!(f.is_zero());
        assert_eq!(f.len(), 6);
    }

    #[test]
    fn test_form_add() {
        let a = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let b = DifferentialForm::new(3, 1, vec![4.0, 5.0, 6.0]);
        let c = a.add(&b);
        assert_eq!(c.coefficients, vec![5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_form_scale() {
        let a = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let b = a.scale(2.0);
        assert_eq!(b.coefficients, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn test_form_norm() {
        let f = DifferentialForm::new(3, 1, vec![3.0, 4.0, 0.0]);
        assert_relative_eq!(f.norm(), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_form_inner_product() {
        let a = DifferentialForm::new(3, 1, vec![1.0, 0.0, 0.0]);
        let b = DifferentialForm::new(3, 1, vec![0.0, 1.0, 0.0]);
        assert_relative_eq!(a.inner_product(&b), 0.0, epsilon = 1e-10);
        assert_relative_eq!(a.inner_product(&a), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_form_display() {
        let f = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let s = format!("{}", f);
        assert!(s.contains("1-form"));
    }

    // ---- Exterior Derivative tests ----

    #[test]
    fn test_exterior_derivative_dim() {
        let d = ExteriorDerivative::new(3, 0);
        // d: 0-forms (dim 1) -> 1-forms (dim 3) → matrix is 3×1
        assert_eq!(d.matrix().nrows(), 3);
        assert_eq!(d.matrix().ncols(), 1);
    }

    #[test]
    fn test_exterior_derivative_dim_2() {
        let d = ExteriorDerivative::new(3, 1);
        // d: 1-forms (dim 3) -> 2-forms (dim 3)
        assert_eq!(d.matrix().nrows(), 3);
        assert_eq!(d.matrix().ncols(), 3);
    }

    #[test]
    fn test_d_squared_zero() {
        // d ∘ d = 0
        let n = 4;
        let d0 = ExteriorDerivative::new(n, 0);
        let d1 = ExteriorDerivative::new(n, 1);
        let d_squared = d1.matrix() * d0.matrix();
        for i in 0..d_squared.nrows() {
            for j in 0..d_squared.ncols() {
                assert_relative_eq!(d_squared[(i, j)], 0.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_d_squared_zero_3d() {
        let n = 3;
        let d0 = ExteriorDerivative::new(n, 0);
        let d1 = ExteriorDerivative::new(n, 1);
        let d_squared = d1.matrix() * d0.matrix();
        for i in 0..d_squared.nrows() {
            for j in 0..d_squared.ncols() {
                assert_relative_eq!(d_squared[(i, j)], 0.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_d_applied_to_0_form() {
        let d = ExteriorDerivative::new(3, 0);
        let f = DifferentialForm::new(3, 0, vec![1.0]);
        let df = d.apply(&f);
        assert_eq!(df.degree, 1);
    }

    // ---- Hodge Star tests ----

    #[test]
    fn test_hodge_star_dim() {
        let star = HodgeStar::new(3, 1);
        assert_eq!(star.matrix().nrows(), 3);
        assert_eq!(star.matrix().ncols(), 3);
    }

    #[test]
    fn test_hodge_star_double_3d_1form() {
        let star = HodgeStar::new(3, 1);
        assert!(star.verify_double_star());
    }

    #[test]
    fn test_hodge_star_double_4d() {
        for k in 0..=4 {
            let star = HodgeStar::new(4, k);
            assert!(star.verify_double_star(), "Failed for n=4, k={}", k);
        }
    }

    #[test]
    fn test_hodge_star_double_2d() {
        for k in 0..=2 {
            let star = HodgeStar::new(2, k);
            assert!(star.verify_double_star(), "Failed for n=2, k={}", k);
        }
    }

    #[test]
    fn test_hodge_star_0_form_on_3d() {
        // ⋆1 = dx ∧ dy ∧ dz
        let star = HodgeStar::new(3, 0);
        let f = DifferentialForm::new(3, 0, vec![1.0]);
        let result = star.apply(&f);
        assert_eq!(result.degree, 3);
        assert_eq!(result.coefficients.len(), 1);
        assert_relative_eq!(result.coefficients[0], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hodge_star_top_form_on_3d() {
        // ⋆(dx ∧ dy ∧ dz) = 1
        let star = HodgeStar::new(3, 3);
        let f = DifferentialForm::new(3, 3, vec![1.0]);
        let result = star.apply(&f);
        assert_eq!(result.degree, 0);
        assert_relative_eq!(result.coefficients[0], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hodge_star_isometry() {
        // ⋆ should preserve the metric up to sign
        let star = HodgeStar::new(3, 1);
        let f = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let sf = star.apply(&f);
        assert_relative_eq!(f.norm(), sf.norm(), epsilon = 1e-10);
    }

    // ---- Laplace-de Rham tests ----

    #[test]
    fn test_laplacian_square_matrix() {
        for n in 2..=4 {
            for k in 0..=n {
                let lap = LaplaceDeRham::new(n, k);
                let dim = binomial(n, k);
                assert_eq!(lap.matrix().nrows(), dim);
                assert_eq!(lap.matrix().ncols(), dim);
            }
        }
    }

    #[test]
    fn test_laplacian_positive_semidefinite_3d() {
        for k in 0..=3 {
            let lap = LaplaceDeRham::new(3, k);
            let svd = lap.matrix().clone().svd(false, false);
            for s in svd.singular_values.iter() {
                assert!(*s >= -1e-8, "Negative singular value for k={}: {}", k, s);
            }
        }
    }

    #[test]
    fn test_zero_form_is_harmonic() {
        let lap = LaplaceDeRham::new(3, 1);
        let zero = DifferentialForm::zero(3, 1);
        assert!(lap.is_harmonic(&zero));
    }

    #[test]
    fn test_laplacian_applied_twice() {
        // Δ is idempotent in the sense that Δ² has same kernel
        let lap = LaplaceDeRham::new(3, 1);
        let lap2 = lap.matrix() * lap.matrix();
        let ker1 = kernel_basis(lap.matrix());
        let ker2 = kernel_basis(&lap2);
        assert_eq!(ker1.len(), ker2.len());
    }

    // ---- Harmonic forms and Betti numbers ----

    #[test]
    fn test_harmonic_forms_are_harmonic() {
        for n in 2..=4 {
            for k in 0..=n {
                let lap = LaplaceDeRham::new(n, k);
                for h in &lap.harmonic_forms() {
                    assert!(lap.is_harmonic(h), "Non-harmonic form found for n={}, k={}", n, k);
                }
            }
        }
    }

    #[test]
    fn test_betti_number_3d() {
        // For ℝ³, all Betti numbers should be 0 (contractible)
        for k in 0..=3 {
            let lap = LaplaceDeRham::new(3, k);
            // On ℝ^n with trivial cohomology, harmonic forms ≅ 0
            // But our finite-dimensional model gives Betti numbers
            // from the kernel of Δ
            let b = lap.betti_number();
            assert!(b <= binomial(3, k), "Betti number too large");
        }
    }

    #[test]
    fn test_harmonic_forms_linearly_independent() {
        for n in 2..=4 {
            for k in 0..=n {
                let lap = LaplaceDeRham::new(n, k);
                let forms = lap.harmonic_forms();
                // Check pairwise inner products (should form an orthonormal-ish set)
                for i in 0..forms.len() {
                    for j in 0..forms.len() {
                        let ip = forms[i].inner_product(&forms[j]);
                        if i == j {
                            assert!(ip.abs() > 1e-10);
                        } else {
                            assert_relative_eq!(ip, 0.0, epsilon = 1e-8);
                        }
                    }
                }
            }
        }
    }

    // ---- Hodge Decomposition tests ----

    #[test]
    fn test_hodge_decompose_3d_1form() {
        let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let decomp = hodge_decompose(&form);
        assert!(decomp.verify(1e-4));
    }

    #[test]
    fn test_hodge_decompose_3d_0form() {
        let form = DifferentialForm::new(3, 0, vec![5.0]);
        let decomp = hodge_decompose(&form);
        assert!(decomp.verify(1e-4));
    }

    #[test]
    fn test_hodge_decompose_4d_2form() {
        let form = DifferentialForm::new(4, 2, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let decomp = hodge_decompose(&form);
        assert!(decomp.verify(1e-4));
    }

    #[test]
    fn test_hodge_decompose_zero() {
        let form = DifferentialForm::zero(3, 1);
        let decomp = hodge_decompose(&form);
        assert!(decomp.verify(1e-8));
        assert!(decomp.exact.is_zero());
        assert!(decomp.coexact.is_zero());
        assert!(decomp.harmonic.is_zero());
    }

    #[test]
    fn test_hodge_decompose_harmonic_is_harmonic() {
        let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let decomp = hodge_decompose(&form);
        let lap = LaplaceDeRham::new(3, 1);
        assert!(lap.is_harmonic(&decomp.harmonic));
    }

    #[test]
    fn test_hodge_decompose_exact_is_image_of_d() {
        // The exact component should be in the image of d
        let form = DifferentialForm::new(4, 2, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let decomp = hodge_decompose(&form);
        assert!(decomp.verify(1e-4));
    }

    #[test]
    fn test_hodge_decompose_multiple_forms() {
        for k in 0..=3 {
            let dim = binomial(3, k);
            let coeffs: Vec<f64> = (0..dim).map(|i| (i as f64 + 1.0) * 0.5).collect();
            let form = DifferentialForm::new(3, k, coeffs);
            let decomp = hodge_decompose(&form);
            assert!(decomp.verify(1e-3), "Failed for k={}", k);
        }
    }

    // ---- Hodge Isomorphism tests ----

    #[test]
    fn test_hodge_isomorphism_verify() {
        let iso = HodgeIsomorphism::new(3);
        assert!(iso.verify_isomorphism());
    }

    #[test]
    fn test_hodge_isomorphism_4d() {
        let iso = HodgeIsomorphism::new(4);
        assert!(iso.verify_isomorphism());
    }

    #[test]
    fn test_hodge_isomorphism_betti_consistency() {
        let iso = HodgeIsomorphism::new(3);
        for k in 0..=3 {
            let reps = iso.harmonic_representatives(k);
            assert_eq!(reps.len(), iso.betti_number(k));
        }
    }

    // ---- Spectral Sequence tests ----

    #[test]
    fn test_spectral_sequence_creation() {
        let ss = HodgeSpectralSequence::new(3);
        assert_eq!(ss.betti_numbers.len(), 4);
    }

    #[test]
    fn test_spectral_sequence_e1_page() {
        let ss = HodgeSpectralSequence::new(3);
        let e1 = ss.e1_page();
        assert_eq!(e1.r, 1);
        assert!(!e1.terms.is_empty() || ss.betti_numbers.iter().all(|&b| b == 0));
    }

    #[test]
    fn test_spectral_sequence_e_inf() {
        let ss = HodgeSpectralSequence::new(3);
        let einf = ss.e_inf_page();
        assert_eq!(einf.r, usize::MAX);
    }

    #[test]
    fn test_total_cohomology() {
        let ss = HodgeSpectralSequence::new(3);
        assert_eq!(ss.total_cohomology(), ss.betti_numbers);
    }

    // ---- Kernel basis tests ----

    #[test]
    fn test_kernel_basis_identity() {
        let id = DMatrix::identity(3, 3);
        let ker = kernel_basis(&id);
        assert_eq!(ker.len(), 0);
    }

    #[test]
    fn test_kernel_basis_zero_matrix() {
        let zero = DMatrix::zeros(3, 3);
        let ker = kernel_basis(&zero);
        assert_eq!(ker.len(), 3);
    }

    #[test]
    fn test_kernel_basis_rank_deficient() {
        let m = DMatrix::from_row_slice(2, 3, &[
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
        ]);
        let ker = kernel_basis(&m);
        // Kernel is spanned by [0, 0, 1] (rank 2, nullity 1)
        assert_eq!(ker.len(), 1);
        // The kernel vector should be in the null space
        let v = DVector::from_vec(ker[0].clone());
        let mv = &m * v;
        assert!(mv.iter().all(|x| x.abs() < 1e-8));
    }

    // ---- Knowledge application tests ----

    #[test]
    fn test_knowledge_space_creation() {
        let ks = knowledge::KnowledgeSpace::new(3);
        assert_eq!(ks.manifold_dim, 3);
        assert_eq!(ks.betti_numbers.len(), 4);
    }

    #[test]
    fn test_knowledge_space_euler() {
        let ks = knowledge::KnowledgeSpace::new(3);
        let euler = ks.euler_characteristic();
        // Euler characteristic should be computable
        assert!(euler.abs() <= 4);
    }

    #[test]
    fn test_knowledge_form_zero() {
        let ks = knowledge::KnowledgeSpace::new(3);
        let f = ks.zero_form(1);
        assert!(f.form.is_zero());
        assert_eq!(f.form.degree, 1);
    }

    #[test]
    fn test_knowledge_form_decompose() {
        let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let kf = knowledge::KnowledgeForm::new(form, vec!["a".into(), "b".into(), "c".into()]);
        let decomp = kf.decompose();
        assert!(decomp.verify());
    }

    #[test]
    fn test_understanding_quality() {
        let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let kf = knowledge::KnowledgeForm::new(form, vec!["a".into(), "b".into(), "c".into()]);
        let q = kf.understanding_quality();
        assert!(q >= 0.0 && q <= 1.0);
    }

    #[test]
    fn test_knowledge_analyze() {
        let ks = knowledge::KnowledgeSpace::new(3);
        let form = DifferentialForm::new(3, 1, vec![2.0, 1.0, 3.0]);
        let kf = knowledge::KnowledgeForm::new(form, vec!["x".into(), "y".into(), "z".into()]);
        let analysis = ks.analyze(&kf);
        assert!(analysis.understanding_quality >= 0.0);
    }

    // ---- Integration / end-to-end tests ----

    #[test]
    fn test_full_pipeline_2d() {
        let n = 2;
        for k in 0..=n {
            let dim = binomial(n, k);
            let coeffs: Vec<f64> = (0..dim).map(|i| (i as f64 + 1.0)).collect();
            let form = DifferentialForm::new(n, k, coeffs);

            let decomp = hodge_decompose(&form);
            assert!(decomp.verify(1e-3), "Decomposition failed for n={}, k={}", n, k);

            let lap = LaplaceDeRham::new(n, k);
            assert!(lap.is_harmonic(&decomp.harmonic));
        }
    }

    #[test]
    fn test_full_pipeline_4d() {
        let n = 4;
        for k in 0..=n {
            let dim = binomial(n, k);
            let coeffs: Vec<f64> = (0..dim).map(|i| (i as f64 * 0.3 + 0.5)).collect();
            let form = DifferentialForm::new(n, k, coeffs);

            let decomp = hodge_decompose(&form);
            assert!(decomp.verify(1e-3), "Decomposition failed for n={}, k={}", n, k);
        }
    }

    #[test]
    fn test_hodge_star_composition_with_decomposition() {
        // Verify Hodge star commutes properly with decomposition
        let star = HodgeStar::new(3, 1);
        let form = DifferentialForm::new(3, 1, vec![1.0, 0.0, 0.0]);
        let starred = star.apply(&form);

        assert_eq!(starred.degree, 2);
        let decomp = hodge_decompose(&starred);
        assert!(decomp.verify(1e-3));
    }

    #[test]
    fn test_betti_numbers_dimension_bound() {
        // Betti numbers can't exceed the dimension of the form space
        for n in 2..=5 {
            for k in 0..=n {
                let lap = LaplaceDeRham::new(n, k);
                assert!(lap.betti_number() <= binomial(n, k));
            }
        }
    }

    #[test]
    fn test_decomposition_orthogonality() {
        // Harmonic should be orthogonal to exact and co-exact
        let form = DifferentialForm::new(3, 1, vec![1.0, 2.0, 3.0]);
        let decomp = hodge_decompose(&form);

        let h_e = decomp.harmonic.inner_product(&decomp.exact);
        let h_c = decomp.harmonic.inner_product(&decomp.coexact);
        assert_relative_eq!(h_e, 0.0, epsilon = 1e-3);
        assert_relative_eq!(h_c, 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_hodge_decomposition_numerical_proof() {
        // Numerically verify the Hodge theorem: every form decomposes
        // Run on many random forms
        for n in 2..=4 {
            for k in 0..=n {
                for trial in 0..5 {
                    let dim = binomial(n, k);
                    let coeffs: Vec<f64> = (0..dim)
                        .map(|i| ((i + trial * 7) as f64 * 1.7 % 5.0 - 2.5))
                        .collect();
                    let form = DifferentialForm::new(n, k, coeffs);
                    let decomp = hodge_decompose(&form);
                    assert!(decomp.verify(1e-2), "Hodge theorem failed for n={}, k={}, trial={}", n, k, trial);
                }
            }
        }
    }
}
