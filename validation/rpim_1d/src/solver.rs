//! Strong-form collocation RPIM solver for 1D steady-state heat conduction.
//!
//! PDE:   -(d/dx) ( k(x) · dT/dx ) = Q(x)   on   x ∈ [x_1, x_N]
//!
//! For *constant* `k`, this reduces to  -k · T''(x) = Q(x).  To keep the code
//! obvious to read, this file handles the constant-k case explicitly. The
//! variable-k generalisation (chain rule gives a term in `k' · T'`) is a short
//! exercise once you've read this.
//!
//! Algorithm (textbook strong-form collocation):
//!
//! 1. Build an `RpimBasis` for the node layout.
//! 2. For every INTERIOR node i, write one equation:
//!        Σ_j [ -k · d²Φ_j/dx²(x_i) ] · T_j = Q(x_i)
//!    (one row of K · T = f).
//! 3. For the two boundary nodes, overwrite the row with the BC:
//!      * Dirichlet T = T_bc            →  row i = e_i,   f_i = T_bc
//!      * Neumann  -k · T'(L) = q_in    →  row i = dΦ/dx, f_i = q_in / (-k) or similar
//! 4. Solve K · T = f with any dense linear solver (we use LU).
//!
//! The `BoundaryCondition::Neumann` convention here is the simple form
//! `dT/dx(x_bc) = g`. If your problem is given as a heat flux `q''` [W/m²],
//! convert first:    g = -q'' / k    (Fourier's law).  An **insulated** wall
//! is just `Neumann(0.0)`.

use nalgebra::{DMatrix, DVector};

use crate::shape::RpimBasis;

/// What to enforce at a boundary node.
#[derive(Debug, Clone, Copy)]
pub enum BoundaryCondition {
    /// Prescribed temperature  T(x_bc) = value
    Dirichlet(f64),
    /// Prescribed temperature gradient  dT/dx(x_bc) = value
    Neumann(f64),
}

/// Everything that defines a steady-state 1D conduction problem.
///
/// `k_const` is used instead of a `k(x)` function so that assembly is a trivial
/// one-liner. Extend to variable `k` by swapping the row expression in
/// [`assemble_and_solve`] for the full chain-rule form.
///
/// Not `Debug`/`Clone` because the `q_source` closure is neither; those
/// derives aren't needed by anything that currently uses this struct.
pub struct Problem {
    pub nodes:    Vec<f64>,
    pub k_const:  f64,
    pub q_source: Box<dyn Fn(f64) -> f64>,
    pub bc_left:  BoundaryCondition,
    pub bc_right: BoundaryCondition,
    /// MQ shape parameter — see docs/theory.md. A reasonable starting value is
    /// 1× the average node spacing.
    pub rbf_c:    f64,
}

impl Problem {
    /// Convenience: build a problem on a uniform grid of `n_nodes` points.
    pub fn uniform(
        n_nodes:  usize,
        x_min:    f64,
        x_max:    f64,
        k_const:  f64,
        q_source: Box<dyn Fn(f64) -> f64>,
        bc_left:  BoundaryCondition,
        bc_right: BoundaryCondition,
    ) -> Self {
        assert!(n_nodes >= 3, "need at least 3 nodes (two BCs + one interior)");
        let dx = (x_max - x_min) / (n_nodes - 1) as f64;
        let nodes = (0..n_nodes).map(|i| x_min + dx * i as f64).collect();
        Self {
            nodes,
            k_const,
            q_source,
            bc_left,
            bc_right,
            // One average node spacing is the engineering default for MQ.
            rbf_c: dx,
        }
    }
}

/// Solved nodal temperatures (one entry per `prob.nodes`).
pub fn assemble_and_solve(prob: &Problem) -> DVector<f64> {
    let n = prob.nodes.len();
    let basis = RpimBasis::new(prob.nodes.clone(), prob.rbf_c);

    let mut k = DMatrix::<f64>::zeros(n, n);
    let mut f = DVector::<f64>::zeros(n);

    // ── Interior collocation (rows 1 .. N-1) ────────────────────────────────
    for i in 1..n - 1 {
        let x_i = prob.nodes[i];
        let s   = basis.evaluate(x_i);
        for j in 0..n {
            k[(i, j)] = -prob.k_const * s.d2phi[j];
        }
        f[i] = (prob.q_source)(x_i);
    }

    // ── Boundary rows: overwrite rows 0 and N-1 ─────────────────────────────
    apply_bc(&mut k, &mut f, 0,     prob.bc_left,  &basis, prob.nodes[0],       prob.k_const);
    apply_bc(&mut k, &mut f, n - 1, prob.bc_right, &basis, prob.nodes[n - 1],   prob.k_const);

    // ── Solve ───────────────────────────────────────────────────────────────
    k.lu().solve(&f)
        .expect("K is singular — check BCs and node layout")
}

/// Overwrite one boundary row of (`k`, `f`) with the chosen BC.
fn apply_bc(
    k:         &mut DMatrix<f64>,
    f:         &mut DVector<f64>,
    row:       usize,
    bc:        BoundaryCondition,
    basis:     &RpimBasis,
    x_boundary:f64,
    k_const:   f64,
) {
    // Clear the row first — RPIM shape functions have global support, so if
    // we leave the pre-existing stiffness contributions in place the BC
    // would compete with them.
    for c in 0..k.ncols() {
        k[(row, c)] = 0.0;
    }

    match bc {
        BoundaryCondition::Dirichlet(t_bc) => {
            // Kronecker-delta property means Φ_j(x_i) = δ_ij, so row = e_i.
            k[(row, row)] = 1.0;
            f[row] = t_bc;
        }
        BoundaryCondition::Neumann(dtdx_bc) => {
            // Σ_j dΦ_j/dx(x_bc) · T_j = dT/dx(x_bc)
            let s = basis.evaluate(x_boundary);
            for (j, &dphi) in s.dphi.iter().enumerate() {
                k[(row, j)] = dphi;
            }
            f[row] = dtdx_bc;
            // `k_const` isn't needed here because we expose the gradient
            // directly (not the heat flux). Kept in the signature for easy
            // extension if a `HeatFlux(q'')` variant is added later.
            let _ = k_const;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_conduction_dirichlet_both_ends() {
        // -k T'' = 0,  T(0) = 100, T(1) = 0  →  T(x) = 100·(1 - x)
        let prob = Problem::uniform(
            11,                                                    // nodes
            0.0, 1.0,
            50.0,                                                  // k
            Box::new(|_x| 0.0),                                    // Q
            BoundaryCondition::Dirichlet(100.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        let t = assemble_and_solve(&prob);
        for (i, &x_i) in prob.nodes.iter().enumerate() {
            let analytic = 100.0 * (1.0 - x_i);
            assert!(
                (t[i] - analytic).abs() < 1e-6,
                "node {i} (x={x_i}): T_rpim = {}, T_analytic = {analytic}",
                t[i],
            );
        }
    }

    #[test]
    fn uniform_heat_generation_dirichlet() {
        // -k T'' = Q₀,  T(0) = 0, T(1) = 0  →  T(x) = Q₀/(2k) · x·(1-x)
        let q0 = 1000.0;
        let k  = 25.0;
        let prob = Problem::uniform(
            15,
            0.0, 1.0,
            k,
            Box::new(move |_x| q0),
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        let t = assemble_and_solve(&prob);
        for (i, &x_i) in prob.nodes.iter().enumerate() {
            let analytic = q0 / (2.0 * k) * x_i * (1.0 - x_i);
            assert!(
                (t[i] - analytic).abs() < 1e-6,
                "node {i} (x={x_i}): T_rpim = {}, T_analytic = {analytic}",
                t[i],
            );
        }
    }
}
