from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import full_proca_solver as f
from relativistic_current_spectrum import initial_parity_control

MODE = ROOT / 'data/raw/linear_handoff_eps1e-3.npz'

# The deterministic handoff is generated on demand so the repository stores
# the calculation, not a derived binary fixture.
if not MODE.exists():
    import subprocess
    subprocess.run(
        [sys.executable, str(ROOT/'src/generate_linear_handoff.py'), '--out', str(MODE)],
        check=True,
        cwd=ROOT,
    )


def test_polarization_helicity_eigenvector():
    for kv in (np.array([1., 2., 3.]), np.array([0., 0., 2.]), np.array([2., -1., .5])):
        kh = kv / np.linalg.norm(kv)
        for s in (+1, -1):
            e = f.polarization(kv, s)
            assert np.allclose(1j*np.cross(kh, e), s*e, rtol=1e-12, atol=1e-12)
            assert abs(np.vdot(kh, e)) < 1e-12


def test_handoff_helicity_label_matches_amplified_sector():
    d = np.load(MODE)
    k = d['kappa']; r = float(d['mass_ratio']); a = float(d['alpha'])
    om2 = k*k + (a*r)**2
    ep = k**3*(np.abs(d['Eplus'])**2 + om2*np.abs(d['Aplus'])**2)
    em = k**3*(np.abs(d['Eminus'])**2 + om2*np.abs(d['Aminus'])**2)
    assert ep.max() > 10*em.max()
    assert np.isclose(float(d['epsilon']), 1.003178166439991e-3, rtol=2e-3)


def test_strict_two_thirds_dealias_excludes_boundary():
    N = 24; L = 6*np.pi
    _, _, dm = f.grids(N, L)
    n = np.fft.fftfreq(N)*N
    idx7 = int(np.where(np.isclose(n, 7))[0][0])
    idx8 = int(np.where(np.isclose(n, 8))[0][0])
    idx0 = int(np.where(np.isclose(n, 0))[0][0])
    assert dm[idx7, idx0, idx0]
    assert not dm[idx8, idx0, idx0]
    assert np.isclose(f.spherical_safe_kmax(N, L), 7/3)


def test_transverse_projection():
    N = 8; K, K2, _ = f.grids(N, 6*np.pi)
    rng = np.random.default_rng(7)
    V = rng.normal(size=(3, N, N, N)) + 1j*rng.normal(size=(3, N, N, N))
    VT = f.project_transverse(V, K, K2)
    dot = sum(K[c]*VT[c] for c in range(3))
    assert np.max(np.abs(dot[K2 > 0])) < 1e-10


def test_gauss_and_independent_lorenz_shadow_smoke(tmp_path):
    _, df = f.run(N=8, L=6*np.pi, dt=.005, tmax=.03, sample=.005,
                  modefile=MODE, out=tmp_path, seed=260823)
    assert np.isfinite(df.select_dtypes(include=[float, int]).to_numpy()).all()
    assert df.gauss_constraint_relative_rms.max() < 1e-12
    assert df.lorenz_shadow_relative_vector_rms.max() < 1e-3
    assert np.all(np.abs(df.chiJ_vector) <= 1 + 1e-10)
    assert (df.longitudinal_field_fraction >= 0).all()


def test_initial_parity_control():
    ctl = initial_parity_control(MODE, [.3, .7], Np=70, Nmu=60)
    assert ctl.chiJ.abs().max() < 1e-8


def test_corrected_headline_tables_exist_and_are_consistent():
    seed = pd.read_csv(ROOT/'data/derived/corrected_seed_robustness.csv')
    assert set(seed.seed) == {260823, 260824, 260825}
    assert (seed.first_main_crossing.between(4.20, 4.30)).all()
    assert (seed.second_main_crossing.between(6.04, 6.08)).all()
    assert (seed.max_positive_chi > 0.40).all()
    assert (seed.endpoint_chi < -0.25).all()

    spatial = pd.read_csv(ROOT/'data/derived/corrected_spatial_convergence.csv')
    assert list(spatial.N) == [24, 32, 40]
    assert (spatial.max_positive_chi > 0.40).all()
    assert (spatial.endpoint_chi < -0.30).all()
    assert spatial.second_main_crossing.max() - spatial.second_main_crossing.min() < 0.04
    assert spatial.common_k_first.notna().all() and spatial.common_k_second.notna().all()

    common = pd.read_csv(ROOT/'data/derived/corrected_common_support_check.csv')
    assert list(common.N) == [25, 32, 40]
    assert np.allclose(common.kinit_effective, 2.4)
    assert common.second_main_crossing.max() - common.second_main_crossing.min() < .05
    assert (common.max_positive_chi > .45).all()

    conv = pd.read_csv(ROOT/'data/derived/corrected_timestep_convergence.csv').iloc[0]
    assert conv.abs_delta_first < 1e-3
    assert conv.abs_delta_second < 1e-3
    assert conv.abs_delta_peak < 5e-4
    assert conv.max_abs_history_delta < 4e-3
    assert conv.fine_max_abs_energy_drift < 5e-5


def test_wigner_conditioning_does_not_create_reversal():
    tab = pd.read_csv(ROOT/'data/derived/corrected_wigner_conditioning.csv')
    assert set(tab.ensemble) == {'conditioned', 'unconditioned'}
    assert (tab.first_main_crossing.notna()).all()
    assert (tab.second_main_crossing.notna()).all()
    assert (tab.max_positive_chi > .4).all()
    assert (tab.endpoint_chi < -.3).all()
    un = tab[tab.ensemble == 'unconditioned'].iloc[0]
    assert np.isclose(un.initial_scale, 1.0)


def test_validation_summary_separates_enforced_and_independent_constraints():
    tab = pd.read_csv(ROOT/'data/derived/corrected_validation_summary.csv').set_index('quantity')
    assert tab.loc['reference_max_gauss_relative', 'value'] < 1e-12
    assert tab.loc['reference_max_lorenz_shadow_relative', 'value'] < 2e-3
    assert tab.loc['reference_max_abs_total_energy_drift', 'value'] < 5e-5
    assert tab.loc['reference_max_longitudinal_fraction', 'value'] > 1e-2


def test_relativistic_current_quadrature_is_realizable():
    tab = pd.read_csv(ROOT/'data/derived/relativistic_current_spectrum.csv')
    assert (tab.P_J_S > 0).all()
    assert (tab.P_J_H.abs() <= tab.P_J_S*(1+1e-8)).all()
    ctl = pd.read_csv(ROOT/'data/derived/relativistic_current_parity_control.csv')
    assert ctl.chiJ.abs().max() < 1e-10


def test_relativistic_current_quadrature_refinement():
    tab = pd.read_csv(ROOT/'data/derived/relativistic_current_quadrature_convergence.csv')
    fine = tab[tab.Np == 320]
    assert set(np.round(fine.K, 2)) == {0.75, 1.50, 2.25}
    assert fine.abs_delta_vs_400.max() < 1.1e-3
    ref = tab[tab.Np == 400]
    assert (ref.P_J_S > 0).all()
    assert (ref.P_J_H.abs() <= ref.P_J_S*(1+1e-10)).all()


def test_dynamic_range_formula_table():
    tab = pd.read_csv(ROOT/'data/derived/corrected_capture_dynamic_range.csv')
    pilot = tab[(tab['case'] == 'corrected pilot endpoint') &
                (tab.UV_margin_times_kpeak == 2.0) &
                (tab.resolved_IR_modes_below_capture == 1)].iloc[0]
    aen = tab[(tab['case'].str.startswith('AEN displayed')) &
              (tab.UV_margin_times_kpeak == 2.0) &
              (tab.resolved_IR_modes_below_capture == 1)].iloc[0]
    assert int(pilot.minimum_N_2over3_dealiased) == 559
    assert int(aen.minimum_N_2over3_dealiased) == 9850
