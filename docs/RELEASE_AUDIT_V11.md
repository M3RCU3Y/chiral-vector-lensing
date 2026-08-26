# Release audit - manuscript v11 / repository v1.1.0

## Executive status

**Analytic source-to-observable chain:** closed at explicitly stated approximation order.

**Production/backreaction pilot:** independently validated and converged across the supplied seed, timestep and spatial tests for the qualitative two-crossing phenomenon.

**Captured-halo absolute lensing forecast:** intentionally open; requires a much larger self-gravitating multiscale simulation.

## Issues discovered and resolved

### 1. Mixed sign conventions

Resolved by fixing one metric/orientation/Fourier/helicity convention and rederiving the physical momentum current. Code tests tie the circular basis to the amplified handoff label.

### 2. Circular Gauss validation

Resolved by changing the paper language: Gauss is algebraically enforced, not an independent validation. Added a non-feedback Lorenz-shadow evolution, total-energy balance, a linear control and convergence tests.

### 3. Exact `N/3` alias boundary

Resolved by strict `|n_i|<N/3` filtering. The old inclusive plane could alias quadratic Chern-Simons products directly back onto retained modes. The independent Lorenz diagnostic exposed the failure even while Gauss remained machine-precision.

### 4. Spatial convergence of the late return

Resolved at pilot level. `24^3`, `32^3`, `40^3` all complete the late negative return. A second `25^3`, `32^3`, `40^3` hierarchy uses exactly common requested `k_init=2.4`, and a common-`k` statistic gives the same qualitative history.

### 5. Finite-volume Wigner normalization

Resolved by deriving the continuum-to-FFT covariance normalization and recording the optional conditioning factor. A same-seed, otherwise identical unconditioned run retains the reversal and late return.

### 6. Relativistic current quadrature

Resolved by replacing a coarse rectangular integral with support-aware Gauss-Legendre quadrature. The release spectrum uses `400x300`; representative `320x240 -> 400x300` refinement changes `chi_J` by at most about `1.03e-3`. Realizability and parity controls are tested.

### 7. Transient lensing versus quasistatic formula

Resolved conceptually. The manuscript introduces unequal-time current spectra and the full-sky transfer structure first. The compact flat-sky/Limber equation is now identified as the quasistatic/equal-time limit rather than advertised as a general transient formula.

### 8. Survey bridge

Added tomographic source kernels and discussion separating gravitational shear from observed ellipticity, including parity-sensitive intrinsic-alignment/systematics considerations.

### 9. Capture-scale overclaim

Resolved with a strict-mask dynamic-range theorem. The paper no longer treats the production box as a capture simulation. With the stated UV margin, one capture-scale interval already needs `N>=559` for the deep pilot and about `N>=9850` for AEN's displayed `537 km/s` escape speed.

### 10. Release package incompleteness

Resolved by adding Makefile targets, requirements, benchmark config, dual licensing, citation/Zenodo metadata, reproducibility docs, generated compact tables and a clean release manifest.

## Acceptance checks for this release

- [x] convention regression tests
- [x] strict de-alias regression test
- [x] deterministic handoff normalization test
- [x] Gauss explicitly labeled enforced
- [x] independent Lorenz-shadow validation
- [x] total-energy drift reported
- [x] three stochastic seeds
- [x] timestep refinement
- [x] full spatial hierarchy through late return
- [x] exactly common-support spatial bridge
- [x] matched linear/no-backreaction control
- [x] matched conditioned/unconditioned Wigner control
- [x] relativistic-current realizability test
- [x] quadrature convergence table
- [x] capture dynamic-range arithmetic test
- [x] manuscript builds with no undefined citations/references
- [x] standalone manuscript PDF visually inspected after rendering (not stored in this repository)

## Deliberately open research tasks

These are not paper defects and should not be disguised as completed calculations:

1. self-gravitating relativistic transport across the full production-to-capture hierarchy;
2. controlled handoff to a nonrelativistic vector halo solver only when `p/m_gamma << 1` is actually satisfied;
3. captured unequal-time `P_J,H^bound(k;eta,eta')`;
4. halo mass/redshift/handedness population model;
5. calibrated survey forecast including intrinsic alignments and instrument/systematic nuisance modeling.

Those form the natural next computational paper rather than unfinished algebra in the present one.
