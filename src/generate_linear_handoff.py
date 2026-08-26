#!/usr/bin/env python3
"""Generate the deterministic AEN vacuum-to-classical handoff.

Repository convention
---------------------
Metric signature: (-,+,+,+).
Conformal field-strength variable: E_i := F_{0i} = A_i' - d_i A_0.
Circular basis: i k_hat x epsilon_s = s epsilon_s.
The translated Chern-Simons sign is chosen so the homogeneous transverse
mode equation used by the nonlinear solver is

    A_s'' + [k^2 + (a r)^2 + s lambda k theta'] A_s = 0.

With theta(0)=+1 this convention calls the dominantly amplified sector s=+1.
Changing the orientation/dual convention exchanges the two circular labels
and flips displayed helicity signs, but no parity-even result.

The stored mode functions include the absolute ma/fa normalization required
for the dimensionless lattice fields A/f_a.  Thus a finite-volume Wigner
sample can be normalized from the continuum covariance rather than by an
arbitrary macroscopic seed.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

GEV_TO_EV = 1.0e9


def alpha_of(tau: float | np.ndarray, hstar: float) -> float | np.ndarray:
    return (1.0 + 0.5*hstar*tau)**2


def hubble_conformal(tau: float | np.ndarray, hstar: float) -> float | np.ndarray:
    a = alpha_of(tau, hstar)
    return hstar/np.sqrt(a)


def mode_frequency_squared(kappa, s, theta_prime, alpha, lam, mass_ratio):
    """Omega_s^2 in the repository circular/sign convention."""
    return kappa*kappa + (alpha*mass_ratio)**2 + s*lam*kappa*theta_prime


def rho_phi_dimensionless(theta, theta_prime, alpha):
    return 0.5*theta_prime**2/alpha**2 + 1.0 - np.cos(theta)


def rho_A_dimensionless(kappa, Aplus, Eplus, Aminus, Eminus, alpha, mass_ratio):
    om_spatial = kappa*kappa + (alpha*mass_ratio)**2
    integrand = (
        np.abs(Eplus)**2 + om_spatial*np.abs(Aplus)**2
        + np.abs(Eminus)**2 + om_spatial*np.abs(Aminus)**2
    )
    # 1/2 field energy and isotropic d^3k/(2pi)^3 = k^2 dk/(2pi^2).
    return 0.5/(2*np.pi**2) * np.trapezoid(kappa*kappa*integrand, kappa) / alpha**4


def rk4_step(tau, dt, theta, pi, Ap, Ep, Am, Em, kappa, lam, r, hstar):
    def deriv(t, th, pp, ap, ep, am, em):
        a = alpha_of(t, hstar)
        hc = hubble_conformal(t, hstar)
        op2 = mode_frequency_squared(kappa, +1, pp, a, lam, r)
        om2 = mode_frequency_squared(kappa, -1, pp, a, lam, r)
        return (
            pp,
            -2*hc*pp - a*a*np.sin(th),
            ep,
            -op2*ap,
            em,
            -om2*am,
        )

    y = (theta, pi, Ap, Ep, Am, Em)
    k1 = deriv(tau, *y)
    y2 = tuple(v + 0.5*dt*k for v, k in zip(y, k1))
    k2 = deriv(tau + 0.5*dt, *y2)
    y3 = tuple(v + 0.5*dt*k for v, k in zip(y, k2))
    k3 = deriv(tau + 0.5*dt, *y3)
    y4 = tuple(v + dt*k for v, k in zip(y, k3))
    k4 = deriv(tau + dt, *y4)
    return tuple(v + dt*(a + 2*b + 2*c + d)/6 for v, a, b, c, d in zip(y, k1, k2, k3, k4))


def generate(
    out: str | Path,
    ma_eV: float = 1e-22,
    fa_GeV: float = 1e17,
    lam: float = 5.0,
    mass_ratio: float = 0.5,
    hstar: float = 7.766231295e-10,
    theta0: float = 1.0,
    pi0: float = 0.0,
    target_epsilon: float = 1e-3,
    kmin: float = 0.02,
    kmax: float = 4.0,
    nk: int = 1000,
    dt: float = 0.005,
    store_interval: float = 0.01,
    taumax: float = 180.0,
):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    k = np.linspace(kmin, kmax, nk)
    a0 = alpha_of(0.0, hstar)
    w0 = np.sqrt(k*k + (a0*mass_ratio)**2)

    # A/f_a covariance in dimensionless coordinates carries ma/fa.
    q = ma_eV/(fa_GeV*GEV_TO_EV)
    A_init = q/np.sqrt(2*w0)
    E_init = -1j*w0*A_init
    Ap = A_init.astype(complex).copy(); Ep = E_init.astype(complex).copy()
    Am = A_init.astype(complex).copy(); Em = E_init.astype(complex).copy()
    Ap0 = Ap.copy(); Ep0 = Ep.copy(); Am0 = Am.copy(); Em0 = Em.copy()
    theta = float(theta0); pi = float(pi0)

    nsteps = int(np.ceil(taumax/dt))
    sample_stride = max(1, int(round(store_interval/dt)))
    crossing = None
    history = []
    for n in range(nsteps + 1):
        tau = n*dt
        if n % sample_stride == 0:
            a = alpha_of(tau, hstar)
            eps = rho_A_dimensionless(k, Ap, Ep, Am, Em, a, mass_ratio) / rho_phi_dimensionless(theta, pi, a)
            history.append((tau, eps, theta, pi, a))
            if eps >= target_epsilon:
                crossing = (tau, eps, theta, pi, a)
                break
        theta, pi, Ap, Ep, Am, Em = rk4_step(
            tau, dt, theta, pi, Ap, Ep, Am, Em, k, lam, mass_ratio, hstar
        )
    if crossing is None:
        raise RuntimeError(f"target epsilon={target_epsilon:g} not reached by tau={taumax:g}")

    tau, eps, theta, pi, a = crossing
    # Occupation computed after removing the absolute ma/fa normalization.
    uAp, uEp = Ap/q, Ep/q
    uAm, uEm = Am/q, Em/q
    w = np.sqrt(k*k + (a*mass_ratio)**2)
    nplus = (np.abs(uEp)**2 + w*w*np.abs(uAp)**2)/(2*w) - 0.5
    nminus = (np.abs(uEm)**2 + w*w*np.abs(uAm)**2)/(2*w) - 0.5

    np.savez_compressed(
        out,
        kappa=k,
        Aplus=Ap, Eplus=Ep, Aminus=Am, Eminus=Em,
        Aplus_initial=Ap0, Eplus_initial=Ep0,
        Aminus_initial=Am0, Eminus_initial=Em0,
        occupation_plus=nplus, occupation_minus=nminus,
        theta=np.array(theta), theta_prime=np.array(pi), alpha=np.array(a),
        tau=np.array(tau), epsilon=np.array(eps), target_epsilon=np.array(target_epsilon),
        mass_ratio=np.array(mass_ratio), lambda_cs=np.array(lam), hstar=np.array(hstar),
        ma_eV=np.array(ma_eV), fa_GeV=np.array(fa_GeV), ma_over_fa=np.array(q),
        metric_signature=np.array("-+++"),
        helicity_convention=np.array("i k_hat x epsilon_s = s epsilon_s"),
        conformal_E_definition=np.array("E_i = F_0i = A_i' - partial_i A_0"),
        cs_sign_convention=np.array("A_s''+[k^2+(a r)^2+s lambda k theta']A_s=0"),
        history=np.asarray(history, float),
    )
    return {
        "path": str(out), "tau": tau, "epsilon": eps, "theta": theta,
        "theta_prime": pi, "alpha": a,
        "occupation_plus_max": float(np.max(nplus)),
        "occupation_minus_max": float(np.max(nminus)),
    }


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="data/raw/linear_handoff_eps1e-3.npz")
    p.add_argument("--ma-eV", type=float, default=1e-22)
    p.add_argument("--fa-GeV", type=float, default=1e17)
    p.add_argument("--lambda-cs", type=float, default=5.0)
    p.add_argument("--mass-ratio", type=float, default=0.5)
    p.add_argument("--target-epsilon", type=float, default=1e-3)
    p.add_argument("--dt", type=float, default=0.005)
    p.add_argument("--nk", type=int, default=1000)
    return p.parse_args()


if __name__ == "__main__":
    a = parse_args()
    info = generate(a.out, ma_eV=a.ma_eV, fa_GeV=a.fa_GeV, lam=a.lambda_cs,
                    mass_ratio=a.mass_ratio, target_epsilon=a.target_epsilon,
                    dt=a.dt, nk=a.nk)
    print("Generated handoff:")
    for k, v in info.items():
        print(f"  {k}: {v}")
