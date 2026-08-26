# Numerical method - v11

## 1. Convention chain

The canonical repository convention is

```text
metric signature: (-,+,+,+)
epsilon^{0123}: +1/sqrt(-g)
Fourier phase: exp(+i k.x)
circular basis: i k_hat x epsilon_s = s epsilon_s
conformal electric variable: E_i = F_0i = A_i' - partial_i A0
B = curl A
```

The orthonormal physical momentum density used for current helicity is

```text
J^hat{i} = T^{hat{0}hat{i}}
         = -a^-4 (E x B)^i - (m_gamma^2/a^2) A0 A_i.
```

Do not replace this by `T_hat0hati` without lowering-index signs.

## 2. Fields and evolution

The nonlinear pilot evolves dimensionless conformal variables

- `theta = phi/f_a`,
- `pi = theta'`,
- spatial Proca potential `A`,
- conformal field-strength variable `E`,
- `A0` reconstructed from Gauss,
- a non-feedback `A0_shadow` evolved from the Lorenz/Proca constraint.

With `alpha=a/a_*`, `r=m_gamma/m_a`, and the release orientation, the implemented first-order system is documented directly in `src/full_proca_solver.py` and derived in the manuscript. The same orientation is used by `generate_linear_handoff.py`, so the handoff helicity labels and nonlinear solver are regression-tested together.

## 3. Vacuum-to-classical initialization

`generate_linear_handoff.py` integrates the transverse circular mode equations from the adiabatic vacuum and stores the absolute `m_a/f_a` normalization. The nonlinear initializer maps the continuum Wigner covariance onto the periodic FFT convention. A finite-volume realization may optionally be conditioned to the target handoff energy; the scale factor is recorded in each output row.

The matched release control keeps the seed, grid, timestep, handoff and UV support fixed and changes only the conditioning switch. The dominant reversal survives essentially unchanged.

## 4. Spatial discretization and strict 2/3 rule

The solver uses a periodic cubic pseudospectral lattice with FFT derivatives and Galerkin de-aliasing of state and nonlinear products.

The retained integer components satisfy

```text
|n_i| < N/3
```

not `<= N/3`. For integer FFT indices this means

```text
|n_i| <= ceil(N/3)-1.
```

The exact boundary is excluded because, when `N` is divisible by three, quadratic products can alias directly back onto the retained boundary. The old inclusive mask could therefore maintain a tiny algebraically enforced Gauss residual while an independent Lorenz diagnostic degraded badly.

For `N=24`, `L=6*pi`, the isotropically safe initialization radius is `7/3`. The common-support convergence bridge therefore uses `N=25,32,40`, for which the requested `k_init=2.4` lies safely inside the retained sphere.

## 5. Time integration

The nonlinear production solver uses explicit midpoint / RK2. The deterministic one-dimensional handoff generator uses RK4.

The release compares `dt=0.02` with `dt=0.01` at `N=24`. The dominant first and second crossing times move by less than `1e-3`; the peak changes by `3.2e-4`; the maximum sampled history difference is about `3.1e-3`.

## 6. Constraint and conservation diagnostics

### Gauss

`A0` is reconstructed algebraically from the nonlinear Proca Gauss equation at every integrator stage. Therefore its residual is an **enforced invariant**, not an independent dynamical validation.

### Lorenz shadow

A shadow `A0_shadow` is evolved using the independent Proca/Lorenz condition and never fed back into the evolution. Agreement with the Gauss-reconstructed `A0` is an independent numerical check. The reference `24^3`, `dt=0.01` run has maximum relative mismatch about `1.06e-3`.

### Energy

The reference maximum absolute total-energy drift is `3.87e-5`.

### Linear control

With nonlinear backreaction disabled, longitudinal excitation is at numerical zero and the Lorenz-shadow residual is about `1e-13`.

## 7. Current helicity

The physical current is Fourier transformed and projected transverse before forming

```text
P_J,S = 1/2 Pi_ij <J_i J_j*>
P_J,H = -i/2 epsilon_ijm k_hat_m <J_i J_j*>
chi_J = P_J,H/P_J,S.
```

The zero mode is excluded. Shell-resolved `chi_J(k,t)` is stored separately for the reference run.

## 8. Relativistic Gaussian current quadrature

The connected Poynting-current covariance is a two-dimensional integral after analytic azimuthal reduction. The current release uses support-aware tensor-product Gauss-Legendre quadrature. For each radial node the angular integration interval is restricted so `q=|K-p|` stays within the tabulated mode-function domain; this avoids artificial discontinuities from zero-padding the interpolant.

Canonical release resolution: `Np=400`, `Nmu=300`.

Representative refinement at `K=0.75,1.5,2.25` is stored in `data/derived/relativistic_current_quadrature_convergence.csv`. The maximum `|delta chi_J|` between `320x240` and `400x300` is about `1.03e-3`. Realizability is an automated test.

## 9. Spatial convergence hierarchy

Two complementary hierarchies are supplied:

1. `24^3,32^3,40^3`: full reference hierarchy, with the `24^3` strict-mask safe sphere ending at `7/3`.
2. `25^3,32^3,40^3`: exactly common requested initialization support `k_init=2.4`.

Both complete the dominant positive excursion and late negative return. A deliberately conservative `k<2.5` current-helicity statistic is also tracked.

## 10. Capture dynamic range

Under the strict mask, a necessary uniform-grid condition for UV production margin `eta_UV`, `n_IR` intervals below capture, production peak `k_pk` and capture momentum `k_cap` is

```text
ceil(N/3)-1 >= eta_UV*n_IR*k_pk/k_cap.
```

Therefore

```text
N_min = 3*ceil(eta_UV*n_IR*k_pk/k_cap)+1.
```

The release table gives `N_min=559` for one capture interval in the deep pilot potential and `N_min=9850` for AEN's displayed `537 km/s` escape speed at the same UV margin. This is a necessary grid criterion, not a sufficient convergence claim.
