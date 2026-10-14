//! Automated validation tests — run with `cargo test --release`.
//!
//! Mirrors the three cases printed by `cargo run`, but as hard assertions so
//! CI can catch regressions in the RPIM shape functions or solver.

use rpim_1d::analytical::{case1_linear, case2_parabolic, case3_insulated_right};
use rpim_1d::solver::{BoundaryCondition, Problem, assemble_and_solve};

const TOL: f64 = 1e-6;
const N_NODES: usize = 21;

fn max_abs_err(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max)
}

#[test]
fn case1_pure_conduction_dirichlet() {
    let (l, k, t_left, t_right) = (1.0, 50.0, 100.0, 0.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(|_| 0.0),
        BoundaryCondition::Dirichlet(t_left),
        BoundaryCondition::Dirichlet(t_right),
    );
    let t_rpim = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case1_linear(x, l, t_left, t_right)).collect();
    let err = max_abs_err(t_rpim.as_slice(), &analytic);
    assert!(err < TOL, "case 1 max error = {err:.3e}");
}

#[test]
fn case2_uniform_heat_generation_dirichlet() {
    let (l, k, q0, t_l, t_r) = (1.0, 25.0, 1000.0, 20.0, 20.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(move |_x| q0),
        BoundaryCondition::Dirichlet(t_l),
        BoundaryCondition::Dirichlet(t_r),
    );
    let t_rpim = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case2_parabolic(x, l, k, q0, t_l, t_r)).collect();
    let err = max_abs_err(t_rpim.as_slice(), &analytic);
    assert!(err < TOL, "case 2 max error = {err:.3e}");
}

#[test]
fn case3_heat_generation_insulated_right() {
    let (l, k, q0, t_left) = (1.0, 25.0, 1000.0, 20.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(move |_x| q0),
        BoundaryCondition::Dirichlet(t_left),
        BoundaryCondition::Neumann(0.0),
    );
    let t_rpim = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case3_insulated_right(x, l, k, q0, t_left)).collect();
    let err = max_abs_err(t_rpim.as_slice(), &analytic);
    assert!(err < TOL, "case 3 max error = {err:.3e}");
}

#[test]
fn refinement_keeps_error_at_machine_precision() {
    // With P2 polynomial augmentation this method reproduces the parabolic
    // steady-state solution of Case 2 EXACTLY (to roundoff) on any grid of
    // ≥3 nodes. We assert that growing the grid never pushes the error above
    // a loose machine-precision band — not strict monotonicity, because at
    // roundoff level the error wobbles non-monotonically between 1e-14 and
    // 1e-11 as the Gram matrix becomes larger and conditioning gets softer.
    let (l, k, q0, t_l, t_r) = (1.0, 25.0, 1000.0, 20.0, 20.0);
    let roundoff_band = 1e-9;
    for &n in &[7usize, 11, 21, 41] {
        let prob = Problem::uniform(
            n, 0.0, l, k,
            Box::new(move |_x| q0),
            BoundaryCondition::Dirichlet(t_l),
            BoundaryCondition::Dirichlet(t_r),
        );
        let t_rpim = assemble_and_solve(&prob);
        let analytic: Vec<f64> = prob.nodes.iter()
            .map(|&x| case2_parabolic(x, l, k, q0, t_l, t_r)).collect();
        let err = max_abs_err(t_rpim.as_slice(), &analytic);
        assert!(
            err < roundoff_band,
            "N={n}: error {err:.3e} exceeded roundoff band {roundoff_band:.0e}",
        );
    }
}
