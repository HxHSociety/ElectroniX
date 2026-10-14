//! Closed-form solutions for the three validation cases.
//!
//! References:
//!   F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine,
//!   "Fundamentals of Heat and Mass Transfer", 7th ed., Wiley 2011.
//!   Chapters 2 (conduction basics) and 3 (one-dimensional, steady-state).
//!
//! Each function below returns the exact `T(x)` for one benchmark problem;
//! the solver's job is to match these to tight tolerance.

/// Case 1 — Dirichlet–Dirichlet, no heat source.
///
///   PDE  : -k · T''(x) = 0
///   BCs  : T(0) = t_left,  T(L) = t_right
///   T(x) : t_left + (t_right - t_left) · x / L
///
/// Incropera §3.2 Eq. 3.6.
pub fn case1_linear(x: f64, l: f64, t_left: f64, t_right: f64) -> f64 {
    t_left + (t_right - t_left) * x / l
}

/// Case 2 — Dirichlet–Dirichlet, uniform volumetric heat generation.
///
///   PDE  : -k · T''(x) = Q0
///   BCs  : T(0) = t_left,  T(L) = t_right
///   T(x) : t_left + (Q0 / (2k)) · x · (L - x)
///                 + (t_right - t_left) · x / L
///
/// Incropera §3.5 Eq. 3.44 (specialised to planar wall with symmetric BCs
/// and extended with linear temperature difference term).
pub fn case2_parabolic(x: f64, l: f64, k: f64, q0: f64, t_left: f64, t_right: f64) -> f64 {
    t_left
        + q0 / (2.0 * k) * x * (l - x)
        + (t_right - t_left) * x / l
}

/// Case 3 — Dirichlet left, insulated right, uniform heat generation.
///
///   PDE  : -k · T''(x) = Q0
///   BCs  : T(0) = t_left,  dT/dx(L) = 0  (perfect insulation at x = L)
///   T(x) : t_left + (Q0 / k) · x · (L - x / 2)
///
/// Incropera §3.5 — plane wall with symmetric half-domain, one wall
/// adiabatic. Also appears in Welty §5.2 Example 5.2.
pub fn case3_insulated_right(x: f64, l: f64, k: f64, q0: f64, t_left: f64) -> f64 {
    t_left + q0 / k * x * (l - 0.5 * x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn case1_boundary_values() {
        assert!((case1_linear(0.0, 1.0, 100.0, 0.0) - 100.0).abs() < 1e-15);
        assert!((case1_linear(1.0, 1.0, 100.0, 0.0) -   0.0).abs() < 1e-15);
    }

    #[test]
    fn case2_matches_known_midpoint_value() {
        // With Q0=1000, k=25, L=1 and both BCs = 0: midpoint should be 5°C.
        let t_mid = case2_parabolic(0.5, 1.0, 25.0, 1000.0, 0.0, 0.0);
        assert!((t_mid - 5.0).abs() < 1e-12);
    }

    #[test]
    fn case3_satisfies_insulation_bc() {
        // Finite-difference derivative at x = L should be ~0.
        let (l, k, q0, t_left) = (1.0, 25.0, 1000.0, 100.0);
        let h = 1e-6;
        let dtdx_at_l = (case3_insulated_right(l, l, k, q0, t_left)
                       - case3_insulated_right(l - h, l, k, q0, t_left)) / h;
        assert!(dtdx_at_l.abs() < 1e-4, "dT/dx(L) = {dtdx_at_l}, want 0");
    }
}
