# Changelog

## v1.1.0 / manuscript v11 - 2026-08-26

### Convention and stress-tensor audit

- fixed one global `(-+++)` metric/orientation/Fourier/helicity convention;
- distinguished the conformal field-strength variable `E_i=F_0i` from the usual orthonormal electric 3-vector;
- rewrote the physical Proca momentum density consistently through the analytic and numerical chain;
- added convention-protecting tests.

### Pseudospectral correction

- replaced the inclusive `|n_i| <= N/3` mask by the strict Orszag rule `|n_i| < N/3`;
- identified exact-boundary aliasing in nonlinear Chern-Simons products as a defect invisible to the algebraically enforced Gauss residual;
- added an independently evolved Lorenz-shadow diagnostic.

### Nonlinear validation

- completed the full late helicity return at `24^3`, `32^3`, and `40^3`;
- added an exactly common-support `25^3`, `32^3`, `40^3` hierarchy;
- retained a conservative common-`k` helicity statistic;
- reran three independent `24^3` seeds;
- completed `dt=0.02 -> 0.01` refinement;
- added a matched no-backreaction control;
- added a matched-seed conditioned/unconditioned Wigner control;
- improved reference total-energy drift to `3.87e-5`.

### Vacuum and Gaussian-current reproducibility

- documented the continuum-to-periodic-box Wigner covariance normalization;
- regenerated the deterministic AEN handoff at `rho_A/rho_phi=1.00318e-3`;
- replaced coarse rectangular current quadrature with support-aware Gauss-Legendre quadrature;
- set the release current spectrum to `Np=400`, `Nmu=300` and added a representative quadrature convergence table;
- added automated realizability and parity controls.

### Lensing and scope

- introduced the unequal-time current covariance required by a transient signal;
- elevated the full-sky transfer structure and demoted the compact flat-sky/Limber equation to its quasistatic/equal-time limit;
- added tomographic-source and intrinsic-alignment/systematics discussion;
- corrected the strict-mask capture lower bound to `N>=559` for the deep pilot and about `N>=9850` for AEN's displayed `537 km/s` escape speed with the stated UV margin;
- explicitly isolates self-gravitating captured-halo `P_J,H` as the remaining multiscale research output.

### Release engineering

- prepared the public companion repository at `M3RCU3Y/chiral-vector-lensing`;
- public GitHub release contains the calculation/data package and omits manuscript source and compiled-paper binaries;
- rebuilt README, Makefile, requirements, licenses, Zenodo/CFF metadata, docs and release audit;
- removed stale provenance-confusing controls and pre-v11 solver backups from the release;
- regenerated corrected figures and verified the manuscript build and test suite.

## v1.0.0 - 2026-08-22

Initial repository-ready internal release. Superseded by v1.1.0 for scientific use.
