//! Validation demo for the 1D RPIM thermal solver.
//!
//! Runs three textbook heat-conduction problems, prints a comparison table
//! (RPIM vs. closed-form), and ends with a pass/fail summary.
//!
//! Usage:
//!     cargo run --release                 # pretty terminal output
//!     cargo run --release -- --csv out    # also dump CSVs to ./out/
//!
//! The `--csv` option writes `case{N}.csv` with columns `x, T_rpim, T_analytic`.
//! Plot them in Python or Excel to see the solution curves.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use rpim_1d::analytical::{case1_linear, case2_parabolic, case3_insulated_right};
use rpim_1d::solver::{BoundaryCondition, Problem, assemble_and_solve};

const TOL_PASS: f64 = 1e-6;   // tolerance below which a case passes
const N_NODES:  usize = 21;   // one node every 5 cm on a 1 m bar

fn main() -> ExitCode {
    let csv_dir = parse_csv_flag();

    println!();
    println!("┌─────────────────────────────────────────────────────────────┐");
    println!("│  rpim_1d — Validation of 1D RPIM Thermal Solver             │");
    println!("│  Reference: Incropera et al., Fundamentals of Heat and      │");
    println!("│             Mass Transfer, 7th ed. (Wiley, 2011) §3.         │");
    println!("└─────────────────────────────────────────────────────────────┘");

    let mut all_passed = true;
    all_passed &= run_case1(&csv_dir);
    all_passed &= run_case2(&csv_dir);
    all_passed &= run_case3(&csv_dir);

    println!();
    println!("─────────────────────────────────────────────────────────────");
    if all_passed {
        println!("  ALL CASES PASSED  (tolerance |err| < {TOL_PASS:.0e})");
        println!("─────────────────────────────────────────────────────────────");
        ExitCode::SUCCESS
    } else {
        println!("  SOME CASES FAILED");
        println!("─────────────────────────────────────────────────────────────");
        ExitCode::FAILURE
    }
}

// ─── CSV option ──────────────────────────────────────────────────────────────

fn parse_csv_flag() -> Option<PathBuf> {
    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        if a == "--csv" {
            let dir = PathBuf::from(args.next().unwrap_or_else(|| "out".into()));
            if let Err(e) = fs::create_dir_all(&dir) {
                eprintln!("warning: cannot create csv dir {}: {}", dir.display(), e);
                return None;
            }
            return Some(dir);
        }
    }
    None
}

fn write_csv(dir: &Path, name: &str, nodes: &[f64], rpim: &[f64], analytic: &[f64]) {
    let p = dir.join(name);
    let mut out = String::from("x,T_rpim,T_analytic,abs_err\n");
    for i in 0..nodes.len() {
        let err = (rpim[i] - analytic[i]).abs();
        out.push_str(&format!("{},{},{},{}\n", nodes[i], rpim[i], analytic[i], err));
    }
    if let Err(e) = fs::write(&p, out) {
        eprintln!("warning: cannot write {}: {}", p.display(), e);
    } else {
        println!("  → wrote {}", p.display());
    }
}

// ─── Pretty-printing ─────────────────────────────────────────────────────────

fn print_header(case_number: usize, title: &str) {
    println!();
    println!("─────────────────────────────────────────────────────────────");
    println!("Case {case_number}  — {title}");
    println!("─────────────────────────────────────────────────────────────");
}

fn print_comparison(nodes: &[f64], rpim: &[f64], analytic: &[f64]) -> f64 {
    println!();
    println!("   {:>8}   {:>13}   {:>13}   {:>10}", "x [m]", "T_RPIM [°C]", "T_analyt [°C]", "|err|");
    let mut max_err = 0.0f64;
    for (i, (&x, (&r, &a))) in nodes.iter().zip(rpim.iter().zip(analytic.iter())).enumerate() {
        let err = (r - a).abs();
        max_err = max_err.max(err);
        // Keep the table short — print every 2nd row.
        if i % 2 == 0 {
            println!("   {x:8.4}   {r:13.6}   {a:13.6}   {err:10.2e}");
        }
    }
    println!();
    let verdict = if max_err < TOL_PASS { "PASS" } else { "FAIL" };
    println!("   max |err| = {max_err:.3e}   →   {verdict}   (tolerance {TOL_PASS:.0e})");
    max_err
}

// ─── Cases ───────────────────────────────────────────────────────────────────

fn run_case1(csv_dir: &Option<PathBuf>) -> bool {
    print_header(1, "Pure conduction, Dirichlet both ends");
    println!("  PDE : -k · d²T/dx² = 0");
    println!("  k   = 50 W/m·K,  L = 1.0 m");
    println!("  BCs : T(0) = 100 °C,   T(L) = 0 °C");
    println!("  Ref : Incropera Eq. 3.6");

    let (l, k, t_left, t_right) = (1.0, 50.0, 100.0, 0.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(|_| 0.0),
        BoundaryCondition::Dirichlet(t_left),
        BoundaryCondition::Dirichlet(t_right),
    );

    let t_rpim    = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case1_linear(x, l, t_left, t_right)).collect();

    let max_err = print_comparison(&prob.nodes, t_rpim.as_slice(), &analytic);
    if let Some(d) = csv_dir { write_csv(d, "case1.csv", &prob.nodes, t_rpim.as_slice(), &analytic); }
    max_err < TOL_PASS
}

fn run_case2(csv_dir: &Option<PathBuf>) -> bool {
    print_header(2, "Uniform volumetric heat generation, Dirichlet both ends");
    println!("  PDE : -k · d²T/dx² = Q₀");
    println!("  k   = 25 W/m·K,  Q₀ = 1000 W/m³,  L = 1.0 m");
    println!("  BCs : T(0) = 20 °C,  T(L) = 20 °C");
    println!("  Ref : Incropera Eq. 3.44 (symmetric Dirichlet form)");

    let (l, k, q0, t_left, t_right) = (1.0, 25.0, 1000.0, 20.0, 20.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(move |_x| q0),
        BoundaryCondition::Dirichlet(t_left),
        BoundaryCondition::Dirichlet(t_right),
    );

    let t_rpim    = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case2_parabolic(x, l, k, q0, t_left, t_right)).collect();

    let max_err = print_comparison(&prob.nodes, t_rpim.as_slice(), &analytic);
    if let Some(d) = csv_dir { write_csv(d, "case2.csv", &prob.nodes, t_rpim.as_slice(), &analytic); }
    max_err < TOL_PASS
}

fn run_case3(csv_dir: &Option<PathBuf>) -> bool {
    print_header(3, "Uniform heat generation, insulated right wall");
    println!("  PDE : -k · d²T/dx² = Q₀");
    println!("  k   = 25 W/m·K,  Q₀ = 1000 W/m³,  L = 1.0 m");
    println!("  BCs : T(0) = 20 °C,  dT/dx(L) = 0 (insulated)");
    println!("  Ref : Incropera §3.5 (plane wall, one adiabatic face)");

    let (l, k, q0, t_left) = (1.0, 25.0, 1000.0, 20.0);
    let prob = Problem::uniform(
        N_NODES, 0.0, l, k,
        Box::new(move |_x| q0),
        BoundaryCondition::Dirichlet(t_left),
        BoundaryCondition::Neumann(0.0),
    );

    let t_rpim    = assemble_and_solve(&prob);
    let analytic: Vec<f64> = prob.nodes.iter()
        .map(|&x| case3_insulated_right(x, l, k, q0, t_left)).collect();

    let max_err = print_comparison(&prob.nodes, t_rpim.as_slice(), &analytic);
    if let Some(d) = csv_dir { write_csv(d, "case3.csv", &prob.nodes, t_rpim.as_slice(), &analytic); }
    max_err < TOL_PASS
}
