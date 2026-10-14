//! Multiquadric Radial Basis Function (MQ-RBF) and its first two derivatives.
//!
//! MQ-RBF:    R(r, c) = sqrt(r² + c²)
//!
//! Here `r = x - x_i` is the SIGNED 1D distance (not |x - x_i|) because we
//! take derivatives with respect to `x`. The RBF value itself is the same
//! for ±r, but its derivatives are not.
//!
//! Hardy's multiquadric is the standard choice for engineering RPIM — it has
//! unbounded support and smooth derivatives, which is exactly what we need
//! for a strong-form collocation solver.
//!
//! References:
//!   G. R. Liu, "Mesh Free Methods: Moving Beyond the Finite Element Method",
//!   CRC Press 2003, §4.4 (MQ-RBF) and §6.3 (RPIM shape functions).

/// Value: R(x; x_i, c) = sqrt((x - x_i)² + c²)
#[inline]
pub fn mq(x: f64, x_i: f64, c: f64) -> f64 {
    let r = x - x_i;
    (r * r + c * c).sqrt()
}

/// First derivative: dR/dx = (x - x_i) / sqrt((x - x_i)² + c²)
#[inline]
pub fn mq_d1(x: f64, x_i: f64, c: f64) -> f64 {
    let r = x - x_i;
    r / (r * r + c * c).sqrt()
}

/// Second derivative: d²R/dx² = c² / ((x - x_i)² + c²)^(3/2)
///
/// Note this is always non-negative (MQ is convex) which is a useful mental
/// check when debugging a bad K matrix.
#[inline]
pub fn mq_d2(x: f64, x_i: f64, c: f64) -> f64 {
    let r = x - x_i;
    let s = r * r + c * c;
    c * c / (s * s.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mq_value_matches_hand_calculation() {
        // At x = x_i the value is just c.
        assert!((mq(3.0, 3.0, 2.0) - 2.0).abs() < 1e-15);
        // Pythagoras: sqrt(3² + 4²) = 5.
        assert!((mq(0.0, 3.0, 4.0) - 5.0).abs() < 1e-15);
    }

    #[test]
    fn mq_derivatives_vs_finite_difference() {
        // Numerically differentiate the RBF and compare — a classic sanity test.
        let c = 0.3;
        let x_i = 1.7;
        let x  = 2.4;
        let h  = 1e-6;

        let d1_num = (mq(x + h, x_i, c) - mq(x - h, x_i, c)) / (2.0 * h);
        assert!((mq_d1(x, x_i, c) - d1_num).abs() < 1e-7);

        let d2_num = (mq(x + h, x_i, c) - 2.0 * mq(x, x_i, c) + mq(x - h, x_i, c)) / (h * h);
        assert!((mq_d2(x, x_i, c) - d2_num).abs() < 1e-4);
    }

    #[test]
    fn mq_second_derivative_is_non_negative() {
        // MQ is convex everywhere; d²R/dx² ≥ 0 for every x.
        for x in [-5.0, -1.0, 0.0, 0.5, 3.0, 10.0] {
            assert!(mq_d2(x, 0.0, 1.0) >= 0.0, "d²R/dx²({x}) was negative");
        }
    }
}
