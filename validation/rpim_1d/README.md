# rpim_1d — Educational 1D RPIM Thermal Solver

A **standalone, dependency-light** Rust project that teaches the Radial Point
Interpolation Method (RPIM) by solving the 1D steady-state heat equation and
comparing against **closed-form textbook answers**.

- **~400 lines of commented code**, no GPU, no external meshes.
- **Only dependency:** [`nalgebra`](https://nalgebra.org/) for linear algebra.
- **Three validation cases** from standard heat-transfer references, each with
  its analytical solution so the error is unambiguous.

## Why does this exist?

ElectroniX's production solver (`rpim_solver` crate) uses RPIM in 3D with GPU
acceleration, local support domains, Anand creep coupling, Darveaux fatigue,
and more. For a new developer that is overwhelming. This project strips
everything away and shows **only** the core of RPIM:

1. Pick a set of scattered nodes.
2. Build *shape functions* at any evaluation point from radial basis functions
   (RBFs) + an augmented polynomial.
3. Discretise the PDE at each interior node (strong-form collocation).
4. Enforce boundary conditions.
5. Solve `K · T = f`.

Once you understand steps 1–5 on 1D paper, the 3D production code reads the same.

## Quick start

```bash
cd validation/rpim_1d
cargo run --release
```

Expected output (truncated):

```
───────────────────────────────────────────────────────────────────────
Case 1  — Pure conduction, Dirichlet both ends     (Incropera Eq. 3.6)
───────────────────────────────────────────────────────────────────────
PDE : -k · d²T/dx² = 0          k = 50 W/m·K,   L = 1.0 m
BCs : T(0) = 100 °C,  T(L) = 0 °C
Analytical : T(x) = 100 · (1 - x/L)

   x         T_RPIM        T_analytic     |err|
  0.000     100.00000      100.00000     0.00e+00
  0.100      90.00000       90.00000     1.42e-14
  …
Max |T_RPIM − T_analytic|  =  2.84e-14  →  PASS  (tol 1e-6)

───────────────────────────────────────────────────────────────────────
Case 2  — Uniform volumetric heat generation       (Incropera Eq. 3.44)
───────────────────────────────────────────────────────────────────────
…
```

Each case ends with **PASS** / **FAIL** against a tight tolerance.

To run the automated tests only:

```bash
cargo test --release
```

To dump CSV so you can plot the solution in Python/Excel:

```bash
cargo run --release -- --csv out.csv
```

## The RPIM method, in two pages

### 1. The approximation

Given `N` nodes with positions `x_1 … x_N` and unknown values `T_1 … T_N`,
RPIM approximates the field anywhere in the domain as

```
T(x) ≈ Σ_i Φ_i(x) · T_i
```

where `Φ_i` are **RPIM shape functions**, built so that `Φ_i(x_j) = δ_ij`
(Kronecker delta). This property is what makes Dirichlet BCs trivial —
just set `T_i` directly.

### 2. Building the shape functions

Each shape function is a linear combination of radial basis functions and
(optionally) low-order polynomials:

```
T(x) = Σ_i a_i · R(||x - x_i||, c)      ← RBF part
      + Σ_k b_k · p_k(x)                  ← polynomial part (we use p_0 = 1, p_1 = x)
```

The RBF this project uses is the **Multiquadric (MQ)**:

```
R(r, c) = sqrt(r² + c²)
```

The constant `c` ("shape parameter") controls how peaked each RBF is;
`c ≈ 1 × average node spacing` is a well-tested default.

Enforcing `T(x_j) = T_j` for every node, plus the standard orthogonality
constraint `Σ_i a_i · p_k(x_i) = 0`, gives the linear system

```
┌ G   P ┐ ┌ a ┐   ┌ T ┐
│       │·│   │ = │   │
└ P^T 0 ┘ └ b ┘   └ 0 ┘
```

where `G[i][j] = R(||x_i - x_j||, c)` is the Gram matrix and `P[i][k] = p_k(x_i)`
is the polynomial basis at the nodes. Call the full augmented matrix `M`.

For any evaluation point `x`, the vector of **shape functions** is

```
Φ(x) = v(x)^T · M^-1       keeping only the first N components
```

with `v(x) = [R(||x - x_1||, c), …, R(||x - x_N||, c), 1, x]`.

Derivatives follow the same formula with `v` replaced by `dv/dx` or `d²v/dx²`:

```
dΦ/dx(x)  = (dv/dx)^T  · M^-1
d²Φ/dx²(x)= (d²v/dx²)^T· M^-1
```

See [`src/shape.rs`](src/shape.rs) for the exact code.

### 3. Discretising the PDE (strong-form collocation)

For 1D steady-state conduction with piecewise-constant thermal conductivity,

```
-k · d²T/dx² = Q(x)
```

We **collocate** the equation at each interior node `x_i`:

```
-k(x_i) · d²T/dx²(x_i) = Q(x_i)
Σ_j [ -k(x_i) · d²Φ_j/dx²(x_i) ] · T_j = Q(x_i)
```

That gives one row of the global linear system `K · T = f` per interior node.

### 4. Boundary conditions

- **Dirichlet** at node `i`: replace row `i` of `K` with the unit vector `e_i`
  and set `f_i = T_prescribed`. (The Kronecker-delta property makes this exact.)
- **Neumann** at node `i`: use the first-derivative shape functions as the row:
  `Σ_j dΦ_j/dx(x_i) · T_j = q_flux / k(x_i)`.

### 5. Solve

`T = K⁻¹ · f` with any linear solver. We use `nalgebra`'s LU decomposition.

## File map

```
validation/rpim_1d/
├── Cargo.toml
├── README.md               ← you are here
├── docs/
│   └── theory.md           ← longer walk-through with references
└── src/
    ├── lib.rs              ← re-exports the public API
    ├── rbf.rs              ← MQ-RBF + derivatives (10 lines of math)
    ├── shape.rs            ← build Φ, dΦ/dx, d²Φ/dx² at an x_eval
    ├── solver.rs           ← assemble K, f, apply BCs, solve
    ├── analytical.rs       ← closed-form reference solutions
    └── main.rs             ← CLI demo that prints the comparison table
tests/
└── validation.rs           ← cargo-test assertions for each case
```

Read them in that order; each file is self-contained and ~80–120 lines.

## Validation cases

| # | Problem | Reference | PDE | Analytical |
|---|---|---|---|---|
| 1 | Fixed-end bar, no source | Incropera 7th ed. §3.2 | `-k·T''  = 0`  | `T(x) = T_0 + (T_L - T_0)·x/L` |
| 2 | Fixed-end bar, uniform heat | Incropera 7th ed. §3.5 | `-k·T''  = Q₀` | `T(x) = T_0 + (Q₀/(2k))·x·(L-x)` |
| 3 | Fixed-left, insulated right | Incropera 7th ed. §3.5 | `-k·T''  = Q₀`, `T'(L) = 0` | `T(x) = T_0 + (Q₀/k)·x·(L - x/2)` |

> **References:**
> F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine,
> *Fundamentals of Heat and Mass Transfer*, 7th ed., Wiley 2011,
> chapters 2–3. (Any heat-transfer textbook has the same formulas.)

## Next steps for curious readers

1. **Change the shape parameter `c`** in `main.rs` — watch what happens for
   very small or very large values (ill-conditioning vs over-smoothing).
2. **Non-uniform node spacing** — RPIM shines here, just change `make_uniform_nodes`
   to any node distribution and the solver still works.
3. **Local support domains** — our solver uses *global* support (every node
   supports every other node) for simplicity. For `N > 100` you want to pick
   only the `k` nearest nodes. See the production `rpim_solver` crate for how.
4. **Transient conduction** — add a time loop with implicit Euler using the same
   shape functions and a mass matrix `M_ij = Φ_i(x_j)` (which equals δ_ij here!).

That last point is the bridge to the production code: once you can solve
steady-state in 1D on paper, the Anand creep solver is just the same ideas in
3D with extra tensor bookkeeping.
