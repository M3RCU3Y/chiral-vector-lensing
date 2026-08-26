PYTHON ?= python
MODE := data/raw/linear_handoff_eps1e-3.npz

.PHONY: test handoff relativistic-current current-convergence postprocess figures smoke release-check reproduce-lite clean

test: handoff
	$(PYTHON) -m pytest -q

handoff:
	$(PYTHON) src/generate_linear_handoff.py --out $(MODE)

relativistic-current:
	$(PYTHON) src/relativistic_current_spectrum.py --modefile $(MODE) --out data/derived/relativistic_current_spectrum.csv --kmin 0.25 --kmax 2.75 --nk 16 --Np 400 --Nmu 300

current-convergence:
	$(PYTHON) scripts/check_current_quadrature.py

postprocess:
	$(PYTHON) scripts/postprocess_results.py

figures: postprocess
	$(PYTHON) scripts/make_figures.py

smoke:
	rm -rf results/smoke
	$(PYTHON) src/full_proca_solver.py --N 8 --dt 0.005 --tmax 0.03 --sample 0.005 --seed 260823 --modefile $(MODE) --out results/smoke

release-check: test postprocess figures

reproduce-lite: handoff relativistic-current current-convergence postprocess figures test

clean:
	rm -rf .pytest_cache src/__pycache__ scripts/__pycache__ tests/__pycache__ results/smoke figures/*.png
