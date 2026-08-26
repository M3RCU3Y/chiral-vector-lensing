# Public repository manifest

This public calculation/data companion contains 46 tracked files.

The manuscript source and compiled PDF are intentionally not version-controlled here. The deterministic handoff is generated from the tracked benchmark/code rather than stored as a binary artifact. Corrected validation histories, derived tables, tests and calculation code are included so the numerical chain can be checked without manuscript files.

All canonical nonlinear reference, timestep, seed, common-support, resolution, Wigner-control, linear-control and shell histories used by the release post-processing are stored directly as CSV files under `results/corrected/`.

## Files

- `.gitattributes`
- `.gitignore`
- `.zenodo.json`
- `CHANGELOG.md`
- `CITATION.cff`
- `LICENSE_CODE`
- `LICENSE_PAPER_DATA`
- `MANIFEST.md`
- `Makefile`
- `README.md`
- `config/aen_halo_benchmark.json`
- `data/derived/corrected_capture_dynamic_range.csv`
- `data/derived/corrected_common_support_check.csv`
- `data/derived/corrected_seed_robustness.csv`
- `data/derived/corrected_spatial_convergence.csv`
- `data/derived/corrected_timestep_convergence.csv`
- `data/derived/corrected_validation_summary.csv`
- `data/derived/corrected_wigner_conditioning.csv`
- `data/derived/relativistic_current_parity_control.csv`
- `data/derived/relativistic_current_quadrature_convergence.csv`
- `data/derived/relativistic_current_spectrum.csv`
- `data/raw/README.md`
- `docs/CLAIMS_MATRIX.md`
- `docs/FIGURE_PROVENANCE.md`
- `docs/NUMERICAL_METHOD.md`
- `docs/RELEASE_AUDIT_V11.md`
- `docs/REPRODUCIBILITY.md`
- `figures/README.md`
- `requirements.txt`
- `results/corrected/linear_control_N24_dt0.02_seed260823.csv`
- `results/corrected/nonlinear_N24_dt0.01_seed260823_full.csv`
- `results/corrected/nonlinear_N24_dt0.02_seed260823_full.csv`
- `results/corrected/nonlinear_N24_dt0.02_seed260823_unconditioned_matched.csv`
- `results/corrected/nonlinear_N24_dt0.02_seed260824_full.csv`
- `results/corrected/nonlinear_N24_dt0.02_seed260825_full.csv`
- `results/corrected/nonlinear_N25_dt0.02_seed260823.csv`
- `results/corrected/nonlinear_N32_dt0.02_seed260823_full.csv`
- `results/corrected/nonlinear_N40_dt0.02_seed260823_full.csv`
- `results/corrected/nonlinear_shells_N24_dt0.01_seed260823_full.csv`
- `scripts/check_current_quadrature.py`
- `scripts/make_figures.py`
- `scripts/postprocess_results.py`
- `src/full_proca_solver.py`
- `src/generate_linear_handoff.py`
- `src/relativistic_current_spectrum.py`
- `tests/test_core.py`
