# Chiral Vector Lensing

**Numerical and reproducibility companion to _Parity-Odd Weak Lensing from Chiral Vector Dark Matter_.**

This repository contains the calculation code, benchmark configuration, corrected numerical histories, derived tables, tests, and figure-generation material used for the manuscript by **Brandon Deonarine**, The University of the West Indies, St. Augustine.

## Physics in one line

The central source variable is the helical part of the **physical transverse momentum-current covariance**,

\[
P_{J,H}(k;\eta,\eta'),
\]

which transfers to vector frame dragging and, in turn, parity-odd weak-lensing shear. The repository implements the relativistic production/current calculations and the corrected nonlinear axion-Proca pilot used to test how that current helicity evolves before the much larger gravitational-capture problem.

## Main numerical result

The corrected nonlinear pilot exhibits

\[
\chi_J<0\;\longrightarrow\;\chi_J>0\;\longrightarrow\;\chi_J<0,
\]

with the full late return surviving:

- three independent `24^3` stochastic seeds;
- timestep refinement `dt=0.02 -> 0.01`;
- `24^3`, `32^3`, `40^3` spatial comparison;
- an exactly common-support `25^3`, `32^3`, `40^3` hierarchy;
- a conservative common-`k` helicity statistic;
- matched no-backreaction and Wigner-conditioning controls;
- strict Orszag `2/3` de-aliasing;
- algebraically enforced Gauss plus an independent Lorenz-shadow diagnostic;
- shell-resolved IR/UV and total-energy checks.

For the reference `24^3`, `dt=0.01`, seed `260823` run,

\[
t_{-\to+}=4.26464,\qquad
\chi_J^{\max}=0.45552,\qquad
t_{+\to-}=6.06734,\qquad
\chi_J(7.5)=-0.37419.
\]

These are **relativistic production/backreaction pilot results**, not a completed prediction for a captured halo population.

## Numerical conventions and correction

The code uses one convention globally,

\[
(-,+,+,+),\qquad
\epsilon^{0123}=+1/\sqrt{-g},\qquad
i\hat{\mathbf k}\times\boldsymbol\epsilon_s=s\boldsymbol\epsilon_s.
\]

The lattice variable is `E_i = F_0i`, and the physical orthonormal Proca momentum density is

\[
J^{\hat i}=T^{\hat0\hat i}
=-a^{-4}(\mathbf E\times\mathbf B)^i
-\frac{m_{\gamma'}^2}{a^2}A_0A_i.
\]

The corrected solver uses the **strict** Orszag `2/3` rule, excluding the exact `|n_i|=N/3` boundary. Keeping that boundary allows quadratic Chern-Simons products to alias directly back onto retained modes while the algebraically reconstructed Gauss residual remains deceptively small. The independent Lorenz-shadow diagnostic detects that failure.

## Repository contents

```text
src/                    production, current-spectrum, nonlinear Proca solvers
scripts/                post-processing, figure and quadrature helpers
tests/                  convention, constraint and convergence tests
config/                 AEN benchmark configuration
data/raw/                deterministic linear handoff covariance
data/derived/            compact tables and quadrature outputs
results/corrected/       canonical corrected nonlinear histories
paper/figures/           generated corrected numerical plots
docs/                    method, claims, provenance and reproducibility notes
```

The manuscript's LaTeX source and compiled PDF are intentionally kept outside this calculation repository. The manuscript itself points to this repository for its public code and numerical data.

## Quick start

Python 3.11+ is recommended.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
make test
```

Useful reproduction targets:

```bash
make handoff                 # deterministic AEN vacuum handoff
make relativistic-current    # 400x300 release current quadrature
make current-convergence     # representative quadrature refinement
make postprocess             # rebuild compact derived tables
make figures                 # regenerate corrected numerical figures
make smoke                   # tiny nonlinear constraint/execution check
make release-check           # tests + tables + figures
```

For the deterministic laptop-level chain including the high-resolution current quadrature:

```bash
make reproduce-lite
```

The full reference nonlinear pilot is intentionally not a routine test:

```bash
python src/full_proca_solver.py \
  --N 24 --dt 0.01 --tmax 7.5 --sample 0.1 --seed 260823 \
  --modefile data/raw/linear_handoff_eps1e-3.npz \
  --save-shells --out results/local-reference
```

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the calculation chain and [`docs/RELEASE_AUDIT_V11.md`](docs/RELEASE_AUDIT_V11.md) for the release-level checks.

## Relativistic Gaussian current calculation

`src/relativistic_current_spectrum.py` evaluates the connected Gaussian Poynting-current covariance

\[
C^J_{ij}(\mathbf K)=\int\!\frac{d^3p}{(2\pi)^3}\,
\epsilon_{iab}\epsilon_{jcd}
\left[C^{EE}_{ac}(\mathbf p)C^{BB}_{bd}(\mathbf q)
+C^{EB}_{ad}(\mathbf p)C^{BE}_{bc}(\mathbf q)\right],
\qquad \mathbf q=\mathbf K-\mathbf p.
\]

The release table uses support-aware Gauss-Legendre quadrature at `Np=400`, `Nmu=300`; representative refinement and the realizability condition `|P_J,H| <= P_J,S` are tested automatically.

## Why gravitational capture is not in this repository

Production and capture occupy an extreme IR/UV hierarchy. With the strict de-alias mask,

\[
N_{\min}=3\left\lceil
\eta_{\rm UV}n_{\rm IR}\frac{k_{\rm pk}}{k_{\rm cap}}
\right\rceil+1.
\]

For the criterion used in the manuscript, the unusually deep pilot potential already requires `N >= 559`; using AEN's displayed `537 km/s` escape speed gives roughly `N >= 9850`. A genuinely self-gravitating production-to-capture calculation is therefore a separate multiscale/HPC problem, not something inferred from the small production box.

## Reproducibility status

The public package distinguishes three levels:

1. **Analytic/numerical checks:** conventions, quadrature, derived tables and test fixtures.
2. **Laptop-scale pilot reproduction:** handoff, current spectrum, post-processing, figures and small nonlinear runs.
3. **Research/HPC follow-up:** self-gravitating production through capture, halo assembly and the final captured unequal-time `P_J,H`.

The third item is the remaining model-dependent physics output. The repository does not manufacture an absolute lensing amplitude without it.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). If using the numerical material, cite the manuscript and this repository.

## License

- Code: MIT, [`LICENSE_CODE`](LICENSE_CODE).
- Original figures and repository-generated data: CC BY 4.0, [`LICENSE_PAPER_DATA`](LICENSE_PAPER_DATA).

Third-party papers are cited in the manuscript and are not redistributed here.
