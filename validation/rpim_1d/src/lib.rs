//! Educational 1D RPIM thermal-conduction solver.
//!
//! Read the README.md next to this file for the full theory walkthrough.
//! Read the modules in this order:
//!
//!   1. [`rbf`]        — the MQ radial basis function and its derivatives
//!   2. [`shape`]      — how we build RPIM shape functions from RBFs
//!   3. [`solver`]     — strong-form collocation assembly of K · T = f
//!   4. [`analytical`] — closed-form answers for the validation cases
//!
//! The `main.rs` binary wires them together into a demo that prints
//! RPIM vs. analytical comparison tables.

pub mod analytical;
pub mod rbf;
pub mod shape;
pub mod solver;

pub use shape::{RpimBasis, Shape};
pub use solver::{BoundaryCondition, Problem, assemble_and_solve};
