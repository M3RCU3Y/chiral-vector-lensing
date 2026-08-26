# Claims matrix - manuscript v11 / repository v1.1.0

This matrix separates theorem-level statements, controlled numerical results, literature inputs and deliberately unresolved physical outputs. It is designed to stop analytic conclusions from quietly absorbing numerical or model-dependent assumptions.

| Claim | Status | Evidence | Scope / limitation |
|---|---|---|---|
| AEN produces unequal transverse circular modes | Literature input + reproduced benchmark | deterministic handoff generator; manuscript production section | linear production stage |
| `m_gamma >> H` does not imply `p/m_gamma << 1` | Analytic | nonrelativistic-handoff proposition | kinematic statement |
| AEN peak modes are initially too fast for immediate ballistic binding in the displayed shallow halo | Analytic / conditional | binding obstruction | necessary condition, not a no-halo theorem |
| Exact transverse relativistic Proca momentum flux is the physical Poynting current in the production sector | Analytic | relativistic stress theorem | `A0=0` transverse stage |
| Statistically isotropic current obeys `|P_J,H| <= P_J,S` | Analytic | current covariance decomposition | no Gaussianity required |
| Parity-symmetric current ensemble has `P_J,H=0` | Analytic | parity proposition | no Gaussianity required |
| Proper-Gaussian relativistic current covariance has the stated convolution form | Analytic | Gaussian current theorem | squeezed/anomalous covariance is additional data |
| Same-helicity relativistic current pairs have a definite parity sign | Analytic | pairwise sign theorem | non-collinear non-null pairs |
| Corrected two-helicity AEN handoff has nonzero realizable current helicity | Numerical quadrature + analytic kernel | `relativistic_current_spectrum.csv` | transverse Gaussian production state |
| Current quadrature is converged at release accuracy | Numerical | `relativistic_current_quadrature_convergence.csv` | representative K values; 320x240 -> 400x300 changes `chi_J` <= 1.03e-3 |
| Vacuum-to-classical lattice covariance has an absolute continuum normalization | Analytic + implementation | Wigner finite-box derivation; handoff metadata | classical-statistical regime |
| Exact finite-volume energy conditioning is not responsible for the reversal | Numerical control | matched conditioned/unconditioned seed-260823 runs | pilot grid |
| Generic nonlinear axion gradients force longitudinal Proca participation | Analytic | nonlinear Gauss proposition | generic inhomogeneous axion field |
| Gauss residual is at machine precision | Enforced invariant | solver diagnostic | **not** an independent validation test |
| Independent Lorenz shadow remains small | Numerical validation | `corrected_validation_summary.csv` | reference max ~1.06e-3 |
| Strict de-aliasing removes the exact-boundary nonlinear alias defect | Numerical/code validation | regression test + Lorenz-shadow comparison | pseudospectral discretization |
| Nonlinear backreaction produces a dominant opposite-sign helicity excursion and late return | Numerical | corrected production runs | pilot production/backreaction problem |
| Dominant reversal/return survives three stochastic seeds | Numerical convergence/control | seed table | `24^3`, `dt=0.02` |
| Dominant reversal/return survives timestep refinement | Numerical convergence | timestep table | `24^3`, `dt=0.02 -> 0.01` |
| Dominant reversal/return survives spatial refinement | Numerical convergence | `24^3,32^3,40^3` table | 24 grid has slightly smaller isotropic UV support |
| Dominant reversal/return survives identical requested UV support | Numerical convergence | `25^3,32^3,40^3` common-support table | all use `k_init=2.4` |
| Conservative common-`k` statistic also reverses and returns | Numerical robustness | spatial/common-support tables | removes strongest UV-edge leverage |
| Matched no-backreaction continuation lacks the large nonlinear positive excursion | Numerical control | linear-control history | same production realization |
| Nonlinear burst creates IR tail and stronger UV cascade | Numerical | shell spectra | resolved pilot band only |
| Resolved vector population remains relativistic during the pilot | Numerical diagnostic | beta-energy history | not a capture simulation |
| Strict uniform-grid production-to-capture lower bound | Analytic | dynamic-range theorem | necessary, not sufficient |
| Leading NR Proca current = convection + spin/magnetization current | Analytic | NR-current theorem | leading gradient order |
| Vector Schrödinger-Poisson conserves particle number and global internal spin | Analytic | conserved-charge proposition | leading scalar SP system |
| Momentum helicity is not conserved by scalar SP gravity | Analytic | helicity-mixing theorem | generic non-collinear scattering |
| Stationary spherical maximally polarized soliton can carry spin/frame dragging but has `P_J,H=0` at displayed order | Analytic | soliton theorems | stationary spherical limit |
| Linear ADM vector constraint maps current helicity to metric-vector helicity | Analytic | vector momentum-constraint theorem | linear metric vector perturbations |
| Unequal-time `P_J,H` maps structurally to full-sky shear `EB` | Analytic structure | full-sky transfer section | Born/linear vector perturbation framework; exact survey kernel conventions stated there |
| Compact flat-sky/Limber formula is a quasistatic/equal-time limit | Analytic approximation statement | lensing section | should not be used for rapidly evolving transient source without unequal-time input |
| Tomographic `C_l,ij^EB` follows after source-distribution weighting | Analytic | survey bridge | observed ellipticity also contains IA/systematics |
| A unique absolute captured-halo `C_l^EB` follows from AEN alone | **False / underdetermined** | identifiability theorem | requires captured unequal-time `P_J,H` + halo population |
| Present pilot proves/disproves gravitational capture | **No** | dynamic-range theorem | required Fourier hierarchy is absent |

## Remaining physical output

The genuinely open quantity is the self-gravitating captured-halo unequal-time spectrum

`P_J,H^bound(k; eta, eta')`

and its halo-population statistics. This requires a much larger relativistic-to-nonrelativistic, self-gravitating initial-value calculation. It is not a missing algebraic lemma and is not assigned a fabricated numerical value in this release.
