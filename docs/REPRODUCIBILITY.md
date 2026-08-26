# Reproducibility guide

This repository is the public calculation/data companion for *Parity-Odd Weak Lensing from Chiral Vector Dark Matter*. The manuscript source and compiled PDF are not stored here.

## Environment

Python 3.11+ is recommended.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Fast verification

```bash
make test
```

The release tests protect the helicity convention, strict de-alias boundary, deterministic handoff normalization, independent Lorenz-shadow diagnostic, seed/timestep/spatial convergence tables, Wigner-conditioning control, relativistic-current realizability and quadrature refinement, and the strict capture-resolution bound.

## Deterministic calculation chain

```bash
make handoff
make relativistic-current
make current-convergence
make postprocess
make figures
make test
```

Equivalent convenience target:

```bash
make reproduce-lite
```

`handoff` regenerates `data/raw/linear_handoff_eps1e-3.npz` from the AEN benchmark. `relativistic-current` evaluates the release `400x300` support-aware Gauss-Legendre current quadrature. `postprocess_results.py` derives the compact validation tables from the checked-in corrected histories rather than hard-coding headline numbers. `make_figures.py` writes regenerated numerical plots to `figures/`.

## Nonlinear reference run

The checked-in reference history can be rerun with

```bash
python src/full_proca_solver.py \
  --N 24 --dt 0.01 --tmax 7.5 --sample 0.1 --seed 260823 \
  --modefile data/raw/linear_handoff_eps1e-3.npz \
  --save-shells --out results/local-reference
```

The larger spatial hierarchy and multi-seed ensemble are research runs rather than CI tasks. Their canonical compact histories are stored directly under `results/corrected/`, and `scripts/postprocess_results.py` reconstructs the release tables from those histories.

## What is not reproduced here

The public repository deliberately stops before the self-gravitating production-to-capture calculation. The manuscript derives why that problem requires a much larger IR/UV hierarchy. The captured unequal-time `P_J,H`, halo population statistics and absolute survey signal therefore remain future multiscale/HPC outputs rather than hidden inputs to this release.
