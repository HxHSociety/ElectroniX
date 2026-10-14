# RPIM Theory — a one-page walk-through

This document complements the inline comments in `src/`. Open it alongside
`src/shape.rs` and `src/solver.rs` for the quickest understanding.

## 1. Scattered-node approximation

Classical finite differences and finite elements both tie their basis to a
**mesh** (cells and element connectivity). Meshless methods — including
**RPIM** — only need a **cloud of points**. There is no element connectivity,
no topology, no cells. This matters for geometry that is hard to mesh
(solder balls on a dense BGA grid) or that changes during the solve
(crack propagation).

Given `N` nodes with positions `x_1 … x_N`, RPIM asks: what is the smoothest
field `T(x)` that interpolates the nodal values `T_1 … T_N`?

## 2. The ansatz

The RPIM approximation at any point `x` is

```
T(x) = Σ_i a_i · R(x; x_i, c)  +  Σ_k b_k · p_k(x)
          \_____________RBFs_____________/   \__polynomials__/
```

- `R(x; x_i, c)` is a radial basis function centred on node `i`. We use the
  **multiquadric** `R(r, c) = sqrt(r² + c²)`.
- `p_k(x)` is a low-order polynomial basis — we use `p_0 = 1`, `p_1 = x`,
  `p_2 = x²`. Including them guarantees the method can reproduce constant,
  linear, **and quadratic** fields exactly — the P2 term matters because
  steady-state conduction with uniform heat generation produces a parabolic
  solution, and anything lower order would fail the patch test (we saw
  ~3 % error in Case 2 before upgrading from P1 to P2).

## 3. The shape parameter `c`

`c` controls how "pointy" the multiquadric is. Classic engineering rule of
thumb:  `c ≈ 1 × average_node_spacing`. Too small and the Gram matrix becomes
near-singular (every RBF looks like a delta). Too large and nearby RBFs
become nearly linearly dependent (also ill-conditioned). The default in
`Problem::uniform` is the exact mean spacing — try changing it in the demo
to feel the sensitivity.

## 4. Augmented linear system

Enforcing `T(x_i) = T_i` for every node AND the constraint
`Σ_i a_i · p_k(x_i) = 0` (needed for solvability) gives the augmented system

```
┌  G    P  ┐  ┌ a ┐    ┌ T ┐
│          │  │   │ =  │   │
└  P^T   0 ┘  └ b ┘    └ 0 ┘
```

Call the full `(N+3) × (N+3)` matrix `M` (three constraint rows because we
chose `{1, x, x²}` as the polynomial basis). Its inverse `M⁻¹` is a *constant*
that depends only on the node layout + the shape parameter + the polynomial
degree. We invert it once in `RpimBasis::new` and reuse it.

## 5. Shape functions and derivatives

Substituting the coefficients back into the ansatz,

```
T(x) = v(x)^T · M⁻¹ · [T; 0]       where v(x) = [R(x;x_1,c), …, R(x;x_N,c), 1, x]^T
     = Φ(x)^T · T                   with Φ(x) = (first N components of M⁻¹^T · v(x))
```

Derivatives of the shape functions follow by differentiating `v` (not `M⁻¹`!):

```
dΦ/dx (x)  = first N components of  M⁻¹^T · dv/dx
d²Φ/dx²(x) = first N components of  M⁻¹^T · d²v/dx²
```

The explicit `v` components are:

| Index   | v(x)                 | dv/dx                | d²v/dx²              |
|---------|----------------------|----------------------|----------------------|
| 1 … N   | `R(x; x_i, c)`       | `(x-x_i)/R`          | `c²/R³`              |
| N+1     | `1`                  | `0`                  | `0`                  |
| N+2     | `x`                  | `1`                  | `0`                  |
| N+3     | `x²`                 | `2x`                 | `2`                  |

## 6. Strong-form collocation

For `-k · T''(x) = Q(x)`, we **collocate** the equation at each interior
node `x_i`:

```
Σ_j [ -k · d²Φ_j/dx²(x_i) ] · T_j = Q(x_i)
```

That's one row of a global linear system `K · T = f` per interior node.
Boundary rows are overwritten:

- **Dirichlet** `T(x_bc) = t`: because of the Kronecker-delta property
  (`Φ_j(x_i) = δ_ij`), the boundary row collapses to a single `1` on the
  diagonal and `t` on the RHS.
- **Neumann** `T'(x_bc) = g`: the boundary row is `Σ_j dΦ_j/dx(x_bc) · T_j = g`.

Dense LU solves `K · T = f`. Done.

## 7. Why this is the production code's teaching sibling

The production `rpim_solver` (3D with Anand creep, GPU PCG, local support
balls, Darveaux fatigue post-processing) uses *exactly* the same shape-
function construction. Local support just means we only pick the `k` nearest
nodes per evaluation point before building `M`, instead of using all `N`.
Everything else in Section 6 generalises component-by-component to 3D
tensors; the solver pipeline (`assemble → BC → solve`) is identical.

## 8. Further reading

- G. R. Liu, *Mesh Free Methods: Moving Beyond the Finite Element Method*,
  CRC Press, 2003. §4 (RBFs), §6 (RPIM).
- G. R. Liu and Y. T. Gu, *An Introduction to Meshfree Methods and Their
  Programming*, Springer, 2005. §6.2–6.3 (RPIM shape functions and
  augmentation).
- F. P. Incropera et al., *Fundamentals of Heat and Mass Transfer*, 7th
  ed., Wiley, 2011. §3 (analytical solutions used for validation).
