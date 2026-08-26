#!/usr/bin/env python3
"""Nonlinear axion-Proca pseudospectral production/backreaction pilot.

Conventions are deliberately explicit because parity signs are load-bearing:

* metric signature (-,+,+,+);
* i k_hat x epsilon_s = s epsilon_s;
* E_i := F_{0i} = A_i' - partial_i A_0 (a conformal field-strength
  variable, the negative of the usual orthonormal electric 3-vector);
* B = curl A;
* the physical orthonormal momentum density is
    J^hat{i} = T^{hat{0}hat{i}}
             = -a^-4 E x B - (r^2/a^2) A_0 A.

The Chern-Simons orientation is the same as ``generate_linear_handoff.py``;
for a homogeneous axion the transverse circular modes obey

    A_s'' + [k^2 + (a r)^2 + s lambda k theta'] A_s = 0.

The solver reconstructs A0 from Gauss at every explicit-midpoint stage.  An
auxiliary ``A0_shadow`` is evolved only with the independent Lorenz/Proca
constraint, A0' = div A - 2 H A0.  It never feeds the evolution and therefore
provides a non-circular constraint diagnostic.

This remains a short relativistic production/backreaction pilot, not a
captured-halo simulation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

DEFAULT_LAMBDA = 5.0
DEFAULT_R = 0.5
DEFAULT_F_OVER_MPL = 1e17 / 2.435e18
DEFAULT_HSTAR = 7.766231295e-10
DEFAULT_ANALYSIS_KMAX = 2.5


def grids(N: int, L: float):
    k = 2*np.pi*np.fft.fftfreq(N, d=L/N)
    K = np.meshgrid(k, k, k, indexing="ij")
    K2 = sum(x*x for x in K)
    n = np.fft.fftfreq(N)*N
    NX, NY, NZ = np.meshgrid(n, n, n, indexing="ij")
    # Strict Orszag 2/3 rule.  The equality mode |n|=N/3 must be excluded
    # when N is divisible by three, otherwise a quadratic product of two
    # boundary modes aliases back onto the retained boundary.
    ncut = int(np.ceil(N/3) - 1)
    dealias = ((np.abs(NX) <= ncut) & (np.abs(NY) <= ncut) & (np.abs(NZ) <= ncut))
    return K, K2, dealias


def spherical_safe_kmax(N: int, L: float) -> float:
    """Largest isotropic spherical k-ball fully inside the componentwise 2/3 mask."""
    return (2*np.pi/L)*int(np.ceil(N/3) - 1)


def filt_real(x, dealias):
    h = np.fft.fftn(x)
    h[~dealias] = 0
    return np.fft.ifftn(h).real


def filt_vec(V, dealias):
    return np.stack([filt_real(V[c], dealias) for c in range(3)])


def grad_hat(h, K):
    return np.stack([np.fft.ifftn(1j*q*h).real for q in K])


def div_vec(V, K, dealias=None):
    Vh = np.stack([np.fft.fftn(V[c]) for c in range(3)])
    if dealias is not None:
        Vh[:, ~dealias] = 0
    return np.fft.ifftn(1j*sum(K[c]*Vh[c] for c in range(3))).real


def project_transverse(Vh, K, K2, zero_mode="zero"):
    out = Vh.copy()
    mask = K2 > 0
    dot = sum(K[c]*out[c] for c in range(3))
    for c in range(3):
        out[c][mask] -= K[c][mask]*dot[mask]/K2[mask]
    if zero_mode == "zero":
        out[:, ~mask] = 0
    elif zero_mode != "preserve":
        raise ValueError("zero_mode must be 'zero' or 'preserve'")
    return out


def polarization(kv, s):
    """Circular polarization satisfying i k_hat x eps_s = s eps_s."""
    kv = np.asarray(kv, float)
    kh = kv/np.linalg.norm(kv)
    ref = np.array([0., 0., 1.]) if abs(kh[2]) <= .9 else np.array([0., 1., 0.])
    e1 = np.cross(ref, kh)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(kh, e1)
    return (e1 + 1j*s*e2)/np.sqrt(2)


def canonical_half_lattice(n):
    for q in n:
        if q != 0:
            return q > 0
    return False


def mode_frequency_squared(k, s, pi, a, lam, mass_ratio):
    """Homogeneous transverse frequency tied to the repository sign convention."""
    return k*k + (a*mass_ratio)**2 + s*lam*k*pi


def physical_momentum_density(E, B, A0, A, a, mass_ratio):
    """T^{hat0 hati}; E here is F_0i, not the usual orthonormal electric vector."""
    poynting_cov = np.moveaxis(
        np.cross(np.moveaxis(E, 0, -1), np.moveaxis(B, 0, -1)), -1, 0
    )
    return -poynting_cov/(a**4) - (mass_ratio**2/(a*a))*A0[None]*A


def operators(theta, A, E, K, K2, dealias, a, lam, mass_ratio):
    th = np.fft.fftn(theta)
    th[~dealias] = 0
    Ah = np.stack([np.fft.fftn(A[c]) for c in range(3)])
    Eh = np.stack([np.fft.fftn(E[c]) for c in range(3)])
    Ah[:, ~dealias] = 0
    Eh[:, ~dealias] = 0

    grad = grad_hat(th, K)
    lap = np.fft.ifftn(-K2*th).real
    KX, KY, KZ = K

    Bh = np.empty_like(Ah)
    Bh[0] = 1j*(KY*Ah[2] - KZ*Ah[1])
    Bh[1] = 1j*(KZ*Ah[0] - KX*Ah[2])
    Bh[2] = 1j*(KX*Ah[1] - KY*Ah[0])
    Bh[:, ~dealias] = 0
    B = np.stack([np.fft.ifftn(Bh[c]).real for c in range(3)])

    divE = np.fft.ifftn(1j*sum(K[c]*Eh[c] for c in range(3))).real
    gdotb = filt_real(np.sum(grad*B, axis=0), dealias)

    # In the chosen dual/orientation convention the Gauss constraint is
    # div E + a^2 r^2 A0 + lambda grad(theta).B = 0.
    A0 = -(divE + lam*gdotb)/((a*mass_ratio)**2)
    A0 = filt_real(A0, dealias)
    A0h = np.fft.fftn(A0); A0h[~dealias] = 0
    gradA0 = grad_hat(A0h, K)

    CBh = np.empty_like(Ah)
    CBh[0] = 1j*(KY*Bh[2] - KZ*Bh[1])
    CBh[1] = 1j*(KZ*Bh[0] - KX*Bh[2])
    CBh[2] = 1j*(KX*Bh[1] - KY*Bh[0])
    CBh[:, ~dealias] = 0
    curlB = np.stack([np.fft.ifftn(CBh[c]).real for c in range(3)])

    return grad, lap, B, curlB, A0, gradA0, Ah, Eh, Bh, divE, gdotb


def rhs(theta, pi, A, E, K, K2, dealias, a, hc, lam, mass_ratio, backreaction=True):
    grad, lap, B, curlB, A0, gradA0, *_ = operators(
        theta, A, E, K, K2, dealias, a, lam, mass_ratio
    )
    EB = filt_real(np.sum(E*B, axis=0), dealias)
    sint = filt_real(np.sin(theta), dealias)
    cross = np.moveaxis(
        np.cross(np.moveaxis(grad, 0, -1), np.moveaxis(E, 0, -1)), -1, 0
    )
    nonlinear = filt_vec(pi[None]*B - cross, dealias)

    dtheta = pi
    dpi = -2*hc*pi + lap - a*a*sint
    if backreaction:
        dpi = dpi + (lam/(a*a))*EB
    dA = E + gradA0
    dE = -curlB - (a*mass_ratio)**2*A - lam*nonlinear
    return dtheta, dpi, dA, dE


def lorenz_shadow_rhs(A, A0_shadow, K, dealias, hc):
    # nabla_mu A^mu=0 -> A0' = div A - 2 Hc A0 for ds^2=a^2(-deta^2+dx^2).
    return filt_real(div_vec(A, K, dealias) - 2*hc*A0_shadow, dealias)


def midpoint(theta, pi, A, E, A0_shadow, dt, K, K2, dealias, a, hc, lam, mass_ratio,
             backreaction=True):
    k1 = rhs(theta, pi, A, E, K, K2, dealias, a, hc, lam, mass_ratio, backreaction)
    s1 = lorenz_shadow_rhs(A, A0_shadow, K, dealias, hc)

    thm = theta + .5*dt*k1[0]
    pim = pi + .5*dt*k1[1]
    Am = A + .5*dt*k1[2]
    Em = E + .5*dt*k1[3]
    A0m = A0_shadow + .5*dt*s1
    k2 = rhs(thm, pim, Am, Em, K, K2, dealias, a, hc, lam, mass_ratio, backreaction)
    s2 = lorenz_shadow_rhs(Am, A0m, K, dealias, hc)

    th = theta + dt*k2[0]
    pp = pi + dt*k2[1]
    AA = A + dt*k2[2]
    EE = E + dt*k2[3]
    A0s = A0_shadow + dt*s2
    return (filt_real(th, dealias), filt_real(pp, dealias),
            filt_vec(AA, dealias), filt_vec(EE, dealias), filt_real(A0s, dealias))


def _complex_interp_log(x, grid, z):
    if x < grid[0] or x > grid[-1]:
        return 0j
    lg = np.log(grid)
    return np.interp(np.log(x), lg, z.real) + 1j*np.interp(np.log(x), lg, z.imag)


def initial_state(N, L, modefile, seed=260823, kinit=2.4, condition_energy=True):
    dat = np.load(modefile)
    kg = dat["kappa"]
    target = float(dat["epsilon"])
    mass_ratio = float(dat["mass_ratio"]) if "mass_ratio" in dat else DEFAULT_R
    lam = float(dat["lambda_cs"]) if "lambda_cs" in dat else DEFAULT_LAMBDA
    modes = {1: (dat["Aplus"], dat["Eplus"]), -1: (dat["Aminus"], dat["Eminus"])}

    Ahat = np.zeros((3, N, N, N), complex)
    Ehat = np.zeros_like(Ahat)
    ksafe = spherical_safe_kmax(N, L)
    kinit_eff = min(float(kinit), float(kg[-1]), ksafe)

    # Continuum -> finite periodic box -> numpy FFT normalization.
    # A(x)=int d^3k/(2pi)^3 A(k)e^{ikx};
    # <|A_box(k)|^2>=V C(k), and numpy ifft divides by N^3.
    fft_covariance_factor = N**3 / (L**1.5)

    for i in range(N):
        ni = i if i < N/2 else i-N
        for j in range(N):
            nj = j if j < N/2 else j-N
            for q in range(N):
                nk = q if q < N/2 else q-N
                nv = (ni, nj, nk)
                if not canonical_half_lattice(nv):
                    continue
                kv = 2*np.pi/L*np.array(nv, float)
                km = np.linalg.norm(kv)
                if km < kg[0] or km > kinit_eff:
                    continue
                par = ((-i) % N, (-j) % N, (-q) % N)
                for s in (1, -1):
                    ss = np.random.SeedSequence([seed, ni+100, nj+100, nk+100, 0 if s == 1 else 1])
                    rg = np.random.default_rng(ss)
                    xi = (rg.normal() + 1j*rg.normal())/np.sqrt(2)
                    ep = polarization(kv, s)
                    Ahat[:, i, j, q] += fft_covariance_factor*xi*_complex_interp_log(km, kg, modes[s][0])*ep
                    Ehat[:, i, j, q] += fft_covariance_factor*xi*_complex_interp_log(km, kg, modes[s][1])*ep
                Ahat[:, par[0], par[1], par[2]] = np.conj(Ahat[:, i, j, q])
                Ehat[:, par[0], par[1], par[2]] = np.conj(Ehat[:, i, j, q])

    K, K2, dealias = grids(N, L)
    Ahat = project_transverse(Ahat, K, K2)
    Ehat = project_transverse(Ehat, K, K2)
    Ahat[:, ~dealias] = 0; Ehat[:, ~dealias] = 0
    A = np.stack([np.fft.ifftn(Ahat[c]).real for c in range(3)])
    E = np.stack([np.fft.ifftn(Ehat[c]).real for c in range(3)])
    theta = np.full((N, N, N), float(dat["theta"]))
    pi = np.full_like(theta, float(dat["theta_prime"]))

    d0 = diagnostics(theta, pi, A, E, K, K2, dealias, 1.0, 0.0, lam, mass_ratio,
                     DEFAULT_F_OVER_MPL, A0_shadow=None)
    pre_eps = d0["epsA"]
    fac = np.sqrt(target/pre_eps) if (condition_energy and pre_eps > 0) else 1.0
    A *= fac; E *= fac
    # Set the independent shadow field only once from the initial Gauss solution.
    *_, A0, _, _, _, _, _, _ = operators(theta, A, E, K, K2, dealias, 1.0, lam, mass_ratio)
    A0_shadow = A0.copy()
    return (theta, pi, A, E, A0_shadow, K, K2, dealias, fac, pre_eps,
            lam, mass_ratio, kinit_eff, ksafe)


def _helicity_from_Jhat(Jh, K, K2, dealias, kmax=None):
    mask = (K2 > 0) & dealias
    if kmax is not None:
        mask &= (K2 <= kmax*kmax)
    km = np.sqrt(K2)
    JT = project_transverse(Jh, K, K2)
    cr = np.moveaxis(
        np.cross(np.moveaxis(JT, 0, -1), np.moveaxis(np.conj(JT), 0, -1)), -1, 0
    )
    H = np.zeros_like(K2, float)
    H[mask] = np.real(-1j*sum(K[c][mask]*cr[c][mask] for c in range(3))/km[mask])
    P = np.sum(np.abs(JT)**2, axis=0)
    pden = P[mask].sum()
    return (float(H[mask].sum()/pden) if pden > 0 else 0.0, H, P, JT, mask)


def diagnostics(theta, pi, A, E, K, K2, dealias, a, hc, lam, mass_ratio, f_over_mpl,
                A0_shadow=None, analysis_kmax=DEFAULT_ANALYSIS_KMAX):
    grad, lap, B, curlB, A0, gradA0, Ah, Eh, Bh, divE, gdotb = operators(
        theta, A, E, K, K2, dealias, a, lam, mass_ratio
    )
    mask = (K2 > 0) & dealias
    km = np.sqrt(K2)

    rho_phi = .5*(pi*pi + np.sum(grad*grad, axis=0))/(a*a) + 1 - np.cos(theta)
    rho_A = .5*(np.sum(E*E, axis=0) + np.sum(B*B, axis=0))/(a**4) \
            + .5*mass_ratio**2*(np.sum(A*A, axis=0) + A0*A0)/(a*a)
    rho = rho_phi + rho_A

    Jv = physical_momentum_density(E, B, A0, A, a, mass_ratio)
    Jh = np.stack([np.fft.fftn(Jv[c]) for c in range(3)])
    chi, H, P, JT, _ = _helicity_from_Jhat(Jh, K, K2, dealias)
    chi_common, *_ = _helicity_from_Jhat(Jh, K, K2, dealias, kmax=analysis_kmax)

    modeE = np.sum(np.abs(Eh)**2 + np.abs(Bh)**2 + (a*mass_ratio)**2*np.abs(Ah)**2, axis=0)
    tot = modeE[mask].sum()
    beta = np.zeros_like(km); beta[mask] = km[mask]/(a*mass_ratio)

    Ath = project_transverse(Ah, K, K2, zero_mode="preserve")
    Eth = project_transverse(Eh, K, K2, zero_mode="preserve")
    longitudinal_power = np.sum((np.abs(Ah-Ath)**2 + np.abs(Eh-Eth)**2)[:, mask])
    nonzero_power = np.sum((np.abs(Ah)**2 + np.abs(Eh)**2)[:, mask])
    longfrac = float(longitudinal_power/nonzero_power) if nonzero_power > 0 else 0.0

    # Algebraically enforced Gauss residual: useful as a code invariant, not an independent test.
    gauss = filt_real(divE + (a*mass_ratio)**2*A0 + lam*gdotb, dealias)
    gauss_scale = np.sqrt(np.mean(divE**2) + np.mean(((a*mass_ratio)**2*A0)**2) + np.mean((lam*gdotb)**2))
    gauss_rel = float(np.sqrt(np.mean(gauss**2))/gauss_scale) if gauss_scale > 0 else 0.0

    # Independent shadow-Lorenz diagnostic.  The shadow does not feed the evolution.
    if A0_shadow is None:
        lorenz_shadow_abs = np.nan
        lorenz_shadow_relvec = np.nan
    else:
        lorenz_shadow_abs = float(np.sqrt(np.mean((A0_shadow-A0)**2)))
        # Normalize the integrated Lorenz-shadow mismatch to the vector-field
        # amplitude, not to A0 itself (which starts exactly at zero and makes a
        # relative A0-only ratio ill-conditioned).
        vecscale = np.sqrt(np.mean(A0*A0) + np.mean(np.sum(A*A, axis=0)))
        lorenz_shadow_relvec = lorenz_shadow_abs/vecscale if vecscale > 0 else 0.0

    rh = np.fft.fftn(rho-rho.mean())
    ph = np.zeros_like(rh)
    ph[mask] = -.5*f_over_mpl**2*a*a*rh[mask]/K2[mask]
    Phi = np.fft.ifftn(ph).real
    vesc = np.sqrt(max(0., -2*Phi.min()))
    vg = beta/np.sqrt(1+beta*beta)

    ksafe = np.max(np.sqrt(K2[(K2 > 0) & dealias &
                               (np.abs(K[0]) <= np.max(np.abs(K[0][dealias])))])) if np.any(mask) else 0.0
    # A practical UV-boundary diagnostic based on the isotropic safe radius of the mask.
    kcomp = min(np.max(np.abs(K[c][dealias])) for c in range(3))
    uv = mask & (km > 0.8*kcomp)
    uvfrac = float(modeE[uv].sum()/tot) if tot > 0 else 0.0

    return dict(
        rho_phi=float(rho_phi.mean()), rho_A=float(rho_A.mean()), rho_total=float(rho.mean()),
        epsA=float(rho_A.mean()/rho_phi.mean()), theta_std=float(theta.std()),
        A0_rms=float(np.sqrt(np.mean(A0*A0))), longitudinal_field_fraction=longfrac,
        gauss_constraint_relative_rms=gauss_rel,
        lorenz_shadow_abs_rms=lorenz_shadow_abs,
        lorenz_shadow_relative_vector_rms=lorenz_shadow_relvec,
        chiJ_vector=chi, chiJ_klt2p5=chi_common,
        beta_energy_mean=float((modeE[mask]*beta[mask]).sum()/tot) if tot > 0 else 0.0,
        frac_beta_lt1=float(modeE[mask & (beta < 1)].sum()/tot) if tot > 0 else 0.0,
        uv_boundary_energy_fraction=uvfrac,
        Phi_rms=float(np.sqrt(np.mean(Phi*Phi))), escape_speed_over_c=float(vesc),
        frac_vgroup_lt_escape=float(modeE[mask & (vg < vesc)].sum()/tot) if (vesc > 0 and tot > 0) else 0.,
        peak_rho=float(rho.max()), rho_std=float(rho.std())
    )


def shell_spectrum(theta, pi, A, E, K, K2, dealias, a, lam, mass_ratio, time):
    grad, lap, B, curlB, A0, gradA0, Ah, Eh, Bh, divE, gdotb = operators(
        theta, A, E, K, K2, dealias, a, lam, mass_ratio
    )
    Jv = physical_momentum_density(E, B, A0, A, a, mass_ratio)
    Jh = np.stack([np.fft.fftn(Jv[c]) for c in range(3)])
    _, H, P, JT, mask = _helicity_from_Jhat(Jh, K, K2, dealias)
    km = np.sqrt(K2)
    dk = min(np.diff(np.unique(np.sort(np.abs(K[0][:,0,0]))))[1:]) if len(np.unique(np.abs(K[0][:,0,0])))>2 else 1.0
    # Fundamental spacing is more robustly 2pi/L, inferred from the first nonzero component.
    comps = np.unique(np.abs(K[0][:,0,0])); comps = comps[comps>0]
    dk = float(comps.min())
    shell = np.rint(km/dk).astype(int)
    mode_energy = np.sum(np.abs(Eh)**2 + np.abs(Bh)**2 + (a*mass_ratio)**2*np.abs(Ah)**2, axis=0)
    rows=[]
    for sh in sorted(np.unique(shell[mask])):
        m = mask & (shell==sh)
        ps=float(P[m].sum()); hh=float(H[m].sum()); ee=float(mode_energy[m].sum())
        rows.append(dict(time=time, shell=int(sh), k_mean=float(km[m].mean()), n_modes=int(m.sum()),
                         P_J_S=ps, P_J_H=hh, chi_J_shell=(hh/ps if ps>0 else 0.0),
                         vector_mode_energy=ee))
    return pd.DataFrame(rows)


def _save_checkpoint(path, th, p, A, E, A0s, time, N, L, dt, seed, fac, pre_eps,
                     lam, mass_ratio, kinit, kinit_eff, ksafe, condition_energy, backreaction):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, theta=th, theta_prime=p, A=A, E=E, A0_shadow=A0s, time=float(time),
        N=int(N), L=float(L), dt=float(dt), seed=int(seed), initial_scale=float(fac),
        precondition_epsA=float(pre_eps), lambda_cs=float(lam), mass_ratio=float(mass_ratio),
        kinit_requested=float(kinit), kinit_effective=float(kinit_eff),
        spherical_safe_kmax=float(ksafe), energy_conditioned=bool(condition_energy),
        backreaction=bool(backreaction))


def _load_checkpoint(path):
    d=np.load(path)
    return d


def run(N=24, L=6*np.pi, dt=.01, tmax=7.5, sample=.25,
        modefile="data/raw/linear_handoff_eps1e-3.npz", out="results", seed=260823,
        hstar=DEFAULT_HSTAR, f_over_mpl=DEFAULT_F_OVER_MPL, kinit=2.4,
        condition_energy=True, backreaction=True, save_shells=False, shell_sample=0.5,
        checkpoint_every=0.0, resume=None):
    modefile = Path(modefile); out = Path(out); out.mkdir(parents=True, exist_ok=True)
    K, K2, dealias = grids(N, L)
    if resume is None:
        (th, p, A, E, A0s, K, K2, dealias, fac, pre_eps, lam, mass_ratio,
         kinit_eff, ksafe) = initial_state(N, L, modefile, seed=seed, kinit=kinit,
                                          condition_energy=condition_energy)
        t0=0.0
    else:
        dcp=_load_checkpoint(resume)
        if int(dcp['N']) != N or abs(float(dcp['L'])-L)>1e-12:
            raise ValueError("checkpoint grid does not match requested N,L")
        if abs(float(dcp['dt'])-dt)>1e-14:
            raise ValueError("checkpoint dt does not match requested dt")
        th=dcp['theta']; p=dcp['theta_prime']; A=dcp['A']; E=dcp['E']; A0s=dcp['A0_shadow']
        t0=float(dcp['time']); seed=int(dcp['seed']); fac=float(dcp['initial_scale'])
        pre_eps=float(dcp['precondition_epsA']); lam=float(dcp['lambda_cs']); mass_ratio=float(dcp['mass_ratio'])
        kinit=float(dcp['kinit_requested']); kinit_eff=float(dcp['kinit_effective']); ksafe=float(dcp['spherical_safe_kmax'])
        condition_energy=bool(dcp['energy_conditioned']); backreaction=bool(dcp['backreaction'])
    rows=[]; shell_rows=[]
    steps=int(round((tmax-t0)/dt)); every=max(1,int(round(sample/dt)))
    shell_every=max(1,int(round(shell_sample/dt)))
    cp_every=max(1,int(round(checkpoint_every/dt))) if checkpoint_every and checkpoint_every>0 else None

    for n in range(steps+1):
        t=t0+n*dt
        a=(1+.5*hstar*t)**2; hc=hstar/np.sqrt(a)
        if n % every == 0 or n == steps:
            d=diagnostics(th,p,A,E,K,K2,dealias,a,hc,lam,mass_ratio,f_over_mpl,A0_shadow=A0s)
            d.update(time=t,N=N,dt=dt,L=L,seed=seed,initial_scale=fac,
                     precondition_epsA=pre_eps,energy_conditioned=bool(condition_energy),
                     lambda_cs=lam,mass_ratio=mass_ratio,kinit_requested=kinit,
                     kinit_effective=kinit_eff,spherical_safe_kmax=ksafe,
                     backreaction=bool(backreaction))
            rows.append(d)
            if save_shells and (n % shell_every == 0 or n == steps):
                shell_rows.append(shell_spectrum(th,p,A,E,K,K2,dealias,a,lam,mass_ratio,t))
            print(f"N={N} t={t:.2f} eps={d['epsA']:.3e} chi={d['chiJ_vector']:+.3f} "
                  f"chi<2.5={d['chiJ_klt2p5']:+.3f} beta={d['beta_energy_mean']:.2f} "
                  f"Long={d['longitudinal_field_fraction']:.2e} Lor={d['lorenz_shadow_relative_vector_rms']:.2e}", flush=True)
        if cp_every is not None and n>0 and n % cp_every == 0:
            tag="nonlinear" if backreaction else "linear_control"
            _save_checkpoint(out/f"{tag}_N{N}_dt{dt:g}_seed{seed}_checkpoint.npz",
                             th,p,A,E,A0s,t,N,L,dt,seed,fac,pre_eps,lam,mass_ratio,
                             kinit,kinit_eff,ksafe,condition_energy,backreaction)
        if n==steps: break
        th,p,A,E,A0s=midpoint(th,p,A,E,A0s,dt,K,K2,dealias,a,hc,lam,mass_ratio,backreaction)

    df=pd.DataFrame(rows)
    # In a resumed segment this is drift relative to the segment start.  The
    # post-processing script also computes a globally anchored version.
    e0=float(df.rho_total.iloc[0])
    df["relative_total_energy_drift"] = df.rho_total/e0 - 1.0
    tag="nonlinear" if backreaction else "linear_control"
    path=out/f"{tag}_N{N}_dt{dt:g}_seed{seed}.csv"; df.to_csv(path,index=False)
    if save_shells:
        sdf=pd.concat(shell_rows,ignore_index=True) if shell_rows else pd.DataFrame()
        sdf.to_csv(out/f"{tag}_shells_N{N}_dt{dt:g}_seed{seed}.csv",index=False)
    return path,df


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N",type=int,default=24); p.add_argument("--L",type=float,default=6*np.pi)
    p.add_argument("--dt",type=float,default=.01); p.add_argument("--tmax",type=float,default=7.5)
    p.add_argument("--sample",type=float,default=.25); p.add_argument("--seed",type=int,default=260823)
    p.add_argument("--kinit",type=float,default=2.4)
    p.add_argument("--modefile",default="data/raw/linear_handoff_eps1e-3.npz")
    p.add_argument("--out",default="results")
    p.add_argument("--unconditioned-wigner",action="store_true")
    p.add_argument("--linear-control",action="store_true")
    p.add_argument("--save-shells",action="store_true")
    p.add_argument("--shell-sample",type=float,default=.5,help="cadence for shell spectra")
    p.add_argument("--checkpoint-every",type=float,default=0.0,help="write restart checkpoint at this cadence")
    p.add_argument("--resume",default=None,help="resume from a solver checkpoint NPZ")
    return p.parse_args()


if __name__=="__main__":
    args=parse_args()
    run(N=args.N,L=args.L,dt=args.dt,tmax=args.tmax,sample=args.sample,
        modefile=args.modefile,out=args.out,seed=args.seed,kinit=args.kinit,
        condition_energy=not args.unconditioned_wigner,
        backreaction=not args.linear_control,save_shells=args.save_shells,
        shell_sample=args.shell_sample,checkpoint_every=args.checkpoint_every,resume=args.resume)
