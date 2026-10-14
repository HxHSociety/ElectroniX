//! RPIM shape functions in 1D with QUADRATIC polynomial augmentation.
//!
//! Given `N` scattered nodes with coordinates `x_1 … x_N`, we build shape
//! functions `Φ_j(x)` such that `Φ_j(x_i) = δ_ij` (Kronecker delta). This
//! property is what makes RPIM so clean for strong-form collocation:
//! Dirichlet boundary conditions become a single-row replacement in the
//! stiffness matrix, with no interpolation error at the node.
//!
//! Why P2 and not P1?
//! ------------------
//! Second-order heat conduction with uniform volumetric heat generation
//! produces a **parabolic** steady-state field (`T(x) ~ x²`). The polynomial
//! augmentation must therefore contain `p_2 = x²` for the method to be
//! *consistent* (reproduce quadratic fields exactly) in strong-form
//! collocation. With P1 only, Cases 2 and 3 in the validation suite fail by
//! a few percent even on fine grids — classic patch-test violation. P2
//! augmentation collapses the error to machine precision.
//!
//! Formulation
//! -----------
//! Approximate a field `u(x)` as a sum of RBFs plus a P2 polynomial:
//!
//! ```text
//! u(x) = Σ_i a_i · R(x; x_i, c) + b_0 · 1 + b_1 · x + b_2 · x²
//! ```
//!
//! Enforcing interpolation at every node + the standard orthogonality
//! constraint (`Σ_i a_i · p_k(x_i) = 0` for each polynomial basis `p_k`)
//! gives the augmented linear system:
//!
//! ```text
//! ┌ G   P  ┐ ┌ a ┐   ┌ u ┐
//! │        │·│   │ = │   │
//! └ P^T 0  ┘ └ b ┘   └ 0 ┘
//!
//! G[i][j] = R(x_i; x_j, c)   (N × N)
//! P[i][k] = p_k(x_i)         (N × 3)   with p_0=1, p_1=x, p_2=x²
//! ```
//!
//! Call the full `(N+3) × (N+3)` matrix `M`. Its inverse maps nodal values
//! `[u; 0]` to the coefficient vector `[a; b]`. Substituting back:
//!
//! ```text
//! u(x) = v(x)^T · M⁻¹ · [u; 0]
//! Φ(x) = v(x)^T · M⁻¹        (we keep only the first N columns)
//! v(x) = [R(x; x_1, c), …, R(x; x_N, c), 1, x, x²]^T
//! ```
//!
//! Derivatives follow immediately — swap `v` for `dv/dx` or `d²v/dx²`:
//!
//! ```text
//! dΦ/dx   = (dv/dx)^T   · M⁻¹
//! d²Φ/dx² = (d²v/dx²)^T · M⁻¹
//! ```
//!
//! Because `M⁻¹` is independent of the evaluation point, we build it ONCE
//! per node layout, store it inside [`RpimBasis`], and reuse it for every
//! `evaluate` call.
//!
//! Reference: Liu & Gu, "An Introduction to Meshfree Methods and Their
//! Programming", Springer 2005, §6.2.

use nalgebra::{DMatrix, DVector};

use crate::rbf::{mq, mq_d1, mq_d2};

/// Reusable RPIM basis for a fixed set of nodes.
///
/// The constructor does the expensive `(N+2)×(N+2)` inversion once.
/// After that, each `evaluate(x_eval)` call is O(N²).
pub struct RpimBasis {
    /// Node coordinates (preserved for RBF value evaluation).
    nodes: Vec<f64>,
    /// MQ-RBF shape parameter.
    c:     f64,
    /// Pre-computed inverse of the augmented Gram matrix.
    /// Rows 0..N correspond to RBF coefficients `a`,
    /// rows N..N+3 correspond to polynomial coefficients `b_0`, `b_1`, `b_2`.
    m_inv: DMatrix<f64>,
}

/// Number of polynomial basis functions used for augmentation: `1, x, x²`.
const N_POLY: usize = 3;

/// Values of `Φ_j(x)`, `dΦ_j/dx`, and `d²Φ_j/dx²` at one evaluation point `x`,
/// one entry per node `j`.
#[derive(Debug, Clone)]
pub struct Shape {
    pub phi:   Vec<f64>,
    pub dphi:  Vec<f64>,
    pub d2phi: Vec<f64>,
}

impl RpimBasis {
    /// Build the basis for `nodes` with MQ shape parameter `c`.
    ///
    /// A good default for `c` is roughly the mean node spacing — pass that in
    /// from the caller so it adapts automatically to coarse vs fine grids.
    pub fn new(nodes: Vec<f64>, c: f64) -> Self {
        let n = nodes.len();
        // Need at least N_POLY nodes so the constraint rows don't drop rank.
        assert!(n >= N_POLY, "need at least {N_POLY} nodes for P2 augmentation");
        assert!(c > 0.0, "MQ shape parameter must be positive");

        let dim = n + N_POLY;
        // Assemble M = [G P; P^T 0] as (N+3) × (N+3).
        let mut m = DMatrix::<f64>::zeros(dim, dim);

        // G[i][j] = R(x_i; x_j, c)
        for i in 0..n {
            for j in 0..n {
                m[(i, j)] = mq(nodes[i], nodes[j], c);
            }
        }

        // P[i] = [1, x_i, x_i²].  Fill both P (columns N..N+3) and P^T (rows).
        for i in 0..n {
            let xi = nodes[i];
            m[(i, n    )] = 1.0;
            m[(i, n + 1)] = xi;
            m[(i, n + 2)] = xi * xi;
            m[(n,     i)] = 1.0;
            m[(n + 1, i)] = xi;
            m[(n + 2, i)] = xi * xi;
        }

        // The lower-right 3×3 block stays zero (constraint rows).
        // Invert — this is the one O((N+3)³) cost of the whole method.
        let m_inv = m
            .try_inverse()
            .expect("augmented Gram matrix is singular — try a different shape parameter c");

        Self { nodes, c, m_inv }
    }

    /// Number of nodes carried by this basis.
    pub fn n(&self) -> usize {
        self.nodes.len()
    }

    /// Evaluate `Φ_j(x_eval)`, `dΦ_j/dx|_{x_eval}`, `d²Φ_j/dx²|_{x_eval}`
    /// for every node index `j`.
    pub fn evaluate(&self, x_eval: f64) -> Shape {
        let n   = self.nodes.len();
        let dim = n + N_POLY;

        // v(x_eval) = [R(x_eval; x_1, c), ..., R(x_eval; x_N, c), 1, x_eval, x_eval²]
        let mut v    = DVector::<f64>::zeros(dim);
        let mut v_d1 = DVector::<f64>::zeros(dim);
        let mut v_d2 = DVector::<f64>::zeros(dim);

        for j in 0..n {
            v[j]    = mq   (x_eval, self.nodes[j], self.c);
            v_d1[j] = mq_d1(x_eval, self.nodes[j], self.c);
            v_d2[j] = mq_d2(x_eval, self.nodes[j], self.c);
        }
        // Polynomial basis contributions:
        //   p_0 = 1   → value 1,      derivatives 0, 0
        //   p_1 = x   → value x_eval, derivatives 1, 0
        //   p_2 = x²  → value x_eval², derivatives 2·x_eval, 2
        v[n]        = 1.0;
        v[n + 1]    = x_eval;
        v[n + 2]    = x_eval * x_eval;
        v_d1[n + 1] = 1.0;
        v_d1[n + 2] = 2.0 * x_eval;
        v_d2[n + 2] = 2.0;

        // Φ(x_eval) is the first N entries of  v^T · M⁻¹ ,
        // which equals (M⁻¹^T · v)[:N] .  Multiplying by m_inv on the LEFT:
        //
        //   row_j(M⁻¹) · v  = Σ_k M⁻¹[j][k] · v[k]
        //
        // Taking only j ∈ 0..N gives Φ_j(x_eval).
        let phi_full   = &self.m_inv * &v;
        let dphi_full  = &self.m_inv * &v_d1;
        let d2phi_full = &self.m_inv * &v_d2;

        Shape {
            phi:   phi_full .as_slice()[..n].to_vec(),
            dphi:  dphi_full.as_slice()[..n].to_vec(),
            d2phi: d2phi_full.as_slice()[..n].to_vec(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(a: f64, b: f64, n: usize) -> Vec<f64> {
        (0..n).map(|i| a + (b - a) * (i as f64) / ((n - 1) as f64)).collect()
    }

    #[test]
    fn kronecker_delta_property() {
        // Φ_j(x_i) should equal δ_ij to high precision.
        let nodes = uniform(0.0, 1.0, 7);
        let basis = RpimBasis::new(nodes.clone(), 0.2);
        for (i, &x_i) in nodes.iter().enumerate() {
            let s = basis.evaluate(x_i);
            for (j, &phi_j) in s.phi.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (phi_j - expected).abs() < 1e-9,
                    "Φ_{j}(x_{i}) = {phi_j}, expected {expected}",
                );
            }
        }
    }

    #[test]
    fn reproduces_constant_field_exactly() {
        // Σ_j Φ_j(x) should equal 1 everywhere (partition of unity),
        // because the polynomial augmentation includes p_0 = 1.
        let basis = RpimBasis::new(uniform(0.0, 1.0, 8), 0.15);
        for frac in [0.0, 0.17, 0.5, 0.83, 1.0] {
            let s   = basis.evaluate(frac);
            let sum: f64 = s.phi.iter().sum();
            assert!((sum - 1.0).abs() < 1e-9, "Σ Φ_j({frac}) = {sum}, want 1");
        }
    }

    #[test]
    fn reproduces_linear_field_exactly() {
        // With p_1 = x in the polynomial basis, RPIM must reproduce u(x) = x
        // exactly: Σ_j Φ_j(x) · x_j = x for every x.
        let nodes = uniform(0.0, 1.0, 6);
        let basis = RpimBasis::new(nodes.clone(), 0.2);
        for frac in [0.1, 0.3, 0.55, 0.9] {
            let s   = basis.evaluate(frac);
            let interp: f64 = s.phi.iter().zip(&nodes).map(|(p, x)| p * x).sum();
            assert!(
                (interp - frac).abs() < 1e-9,
                "Σ Φ_j({frac}) · x_j = {interp}, want {frac}",
            );
        }
    }

    #[test]
    fn reproduces_quadratic_field_exactly() {
        // With p_2 = x² in the polynomial basis, RPIM must reproduce u(x) = x²
        // and give the correct second derivative (2) everywhere.
        let nodes = uniform(0.0, 1.0, 8);
        let basis = RpimBasis::new(nodes.clone(), 0.15);
        let u: Vec<f64> = nodes.iter().map(|x| x * x).collect();
        for frac in [0.0, 0.1, 0.37, 0.5, 0.83, 1.0] {
            let s = basis.evaluate(frac);
            let u_val:    f64 = s.phi  .iter().zip(&u).map(|(p, v)| p * v).sum();
            let u_second: f64 = s.d2phi.iter().zip(&u).map(|(p, v)| p * v).sum();
            assert!(
                (u_val - frac * frac).abs() < 1e-9,
                "u({frac}) = {u_val}, want {}",
                frac * frac,
            );
            assert!(
                (u_second - 2.0).abs() < 1e-6,
                "d²u/dx²({frac}) = {u_second}, want 2",
            );
        }
    }

    #[test]
    fn derivatives_of_linear_field_are_correct() {
        // For u(x) = x,  du/dx = 1 and d²u/dx² = 0, exactly, at every x.
        let nodes = uniform(0.0, 1.0, 6);
        let basis = RpimBasis::new(nodes.clone(), 0.2);
        for frac in [0.1, 0.4, 0.7] {
            let s = basis.evaluate(frac);
            let du:  f64 = s.dphi .iter().zip(&nodes).map(|(p, x)| p * x).sum();
            let d2u: f64 = s.d2phi.iter().zip(&nodes).map(|(p, x)| p * x).sum();
            assert!((du  - 1.0).abs() < 1e-7, "du/dx  at {frac} = {du}");
            assert!( d2u.abs()        < 1e-7, "d²u/dx² at {frac} = {d2u}");
        }
    }
}
