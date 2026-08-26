#!/usr/bin/env python3
"""Relativistic Gaussian Poynting-current covariance for the AEN handoff.

Uses the exact connected equal-time Gaussian formula

 C^J_ij(K)=∫ d^3p/(2π)^3 eps_iab eps_jcd
 [C^EE_ac(p) C^BB_bd(q) + C^EB_ad(p) C^BE_bc(q)], q=K-p,

for the transverse linear state.  The computation is basis-independent once
written in circular projectors.  Absolute normalization cancels from chi_J,
but the stored handoff retains it.
"""
from __future__ import annotations
import argparse
import gc
from pathlib import Path
import numpy as np
import pandas as pd

EPS3 = np.zeros((3,3,3), int)
EPS3[0,1,2]=EPS3[1,2,0]=EPS3[2,0,1]=1
EPS3[0,2,1]=EPS3[2,1,0]=EPS3[1,0,2]=-1


def projector(khat, s):
    khat=np.asarray(khat,float)
    Pi=np.eye(3)-np.outer(khat,khat)
    anti=np.einsum('ijm,m->ij',EPS3,khat)
    return 0.5*(Pi-1j*s*anti)


def interp_complex(x, grid, z):
    if x < grid[0] or x > grid[-1]: return 0j
    return np.interp(x,grid,z.real)+1j*np.interp(x,grid,z.imag)


def covariances(kvec, grid, modes, initial=False):
    k=float(np.linalg.norm(kvec))
    if k==0 or k<grid[0] or k>grid[-1]:
        z=np.zeros((3,3),complex); return z,z,z,z
    kh=kvec/k
    CEE=np.zeros((3,3),complex); CBB=np.zeros_like(CEE)
    CEB=np.zeros_like(CEE); CBE=np.zeros_like(CEE)
    for s in (+1,-1):
        Aarr,Earr=modes[s]
        A=interp_complex(k,grid,Aarr); E=interp_complex(k,grid,Earr)
        P=projector(kh,s)
        Bs=s*k*A
        CEE += (abs(E)**2)*P
        CBB += (abs(Bs)**2)*P
        CEB += (E*np.conj(Bs))*P
        CBE += (Bs*np.conj(E))*P
    return CEE,CBB,CEB,CBE


def _projector_batch(kvec, s):
    """Circular projector for an array (...,3) of nonzero vectors."""
    kvec=np.asarray(kvec,float)
    km=np.linalg.norm(kvec,axis=-1)
    kh=np.zeros_like(kvec)
    good=km>0
    kh[good]=kvec[good]/km[good,None]
    eye=np.eye(3)
    Pi=eye-np.einsum('...i,...j->...ij',kh,kh)
    anti=np.einsum('ijm,...m->...ij',EPS3,kh)
    P=0.5*(Pi-1j*s*anti)
    P[~good]=0
    return P,km


def _interp_complex_array(x, grid, z):
    x=np.asarray(x,float)
    out=np.interp(x,grid,z.real,left=0.0,right=0.0)+1j*np.interp(x,grid,z.imag,left=0.0,right=0.0)
    out[(x<grid[0])|(x>grid[-1])]=0.0
    return out


def _covariances_batch(kvec, grid, modes):
    """Return CEE,CBB,CEB,CBE for kvec with shape (...,3)."""
    shape=kvec.shape[:-1]+(3,3)
    CEE=np.zeros(shape,complex); CBB=np.zeros(shape,complex)
    CEB=np.zeros(shape,complex); CBE=np.zeros(shape,complex)
    for s in (+1,-1):
        P,km=_projector_batch(kvec,s)
        Aarr,Earr=modes[s]
        A=_interp_complex_array(km,grid,Aarr)
        E=_interp_complex_array(km,grid,Earr)
        Bs=s*km*A
        CEE += (np.abs(E)**2)[...,None,None]*P
        CBB += (np.abs(Bs)**2)[...,None,None]*P
        CEB += (E*np.conj(Bs))[...,None,None]*P
        CBE += (Bs*np.conj(E))[...,None,None]*P
    return CEE,CBB,CEB,CBE


def current_covariance(modefile, K, Np=180, Nmu=140, use_initial=False):
    """Connected Gaussian current covariance for ``K`` along z.

    The azimuthal integral is analytic for the transverse trace and the
    z-directed antisymmetric contraction, leaving a two-dimensional integral
    over ``p`` and ``mu=cos(theta)``.  We use tensor-product Gauss-Legendre
    quadrature, but crucially restrict the angular interval at every radial
    node to the domain for which ``q=|K-p|`` lies inside the tabulated mode
    support.  This removes the artificial endpoint discontinuities produced by
    padding the interpolant with zeros and converges substantially faster than
    a rectangular trapezoidal grid.
    """
    d=np.load(modefile)
    grid=np.asarray(d['kappa'],float)
    if use_initial and 'Aplus_initial' in d:
        modes={+1:(d['Aplus_initial'],d['Eplus_initial']),-1:(d['Aminus_initial'],d['Eminus_initial'])}
    else:
        modes={+1:(d['Aplus'],d['Eplus']),-1:(d['Aminus'],d['Eminus'])}

    K=float(K)
    if K <= 0:
        raise ValueError('K must be positive')
    pmin,pmax=float(grid[0]),float(grid[-1])
    xp,wp0=np.polynomial.legendre.leggauss(int(Np))
    xm,wm0=np.polynomial.legendre.leggauss(int(Nmu))
    p=0.5*(pmax-pmin)*xp+0.5*(pmax+pmin)
    wp=0.5*(pmax-pmin)*wp0

    # q^2=K^2+p^2-2 K p mu.  Restrict mu so pmin <= q <= pmax.
    lo=(K*K+p*p-pmax*pmax)/(2*K*p)
    hi=(K*K+p*p-pmin*pmin)/(2*K*p)
    lo=np.maximum(-1.0,lo); hi=np.minimum(1.0,hi)
    width=np.maximum(0.0,hi-lo)
    center=0.5*(hi+lo); half=0.5*width
    mu=center[:,None]+half[:,None]*xm[None,:]
    wmu=half[:,None]*wm0[None,:]

    pp=p[:,None]
    st=np.sqrt(np.maximum(0.,1-mu*mu))
    pvec=np.stack([pp*st,np.zeros_like(mu),pp*mu],axis=-1)
    qvec=np.zeros_like(pvec); qvec[...,2]=K; qvec-=pvec

    EEp,BBp,EBp,BEp=_covariances_batch(pvec,grid,modes)
    EEq,BBq,EBq,BEq=_covariances_batch(qvec,grid,modes)
    term=np.einsum('iab,jcd,...ac,...bd->...ij',EPS3,EPS3,EEp,BBq,optimize=True)
    term+=np.einsum('iab,jcd,...ad,...bc->...ij',EPS3,EPS3,EBp,BEq,optimize=True)
    weight=((2*np.pi)/(2*np.pi)**3)*(pp*pp)*wp[:,None]*wmu
    C=np.einsum('...,...ij->ij',weight,term,optimize=True)
    C=0.5*(C+C.conj().T)
    kh=np.array([0.,0.,1.])
    Pi=np.eye(3)-np.outer(kh,kh)
    PS=0.5*np.real(np.einsum('ij,ij->',Pi,C))
    PH=np.real(-0.5j*np.einsum('ijm,m,ij->',EPS3,kh,C))
    return PS,PH,C

def spectrum(modefile, K_values, Np=180, Nmu=140, use_initial=False):
    rows=[]
    for K in K_values:
        PS,PH,_=current_covariance(modefile,float(K),Np=Np,Nmu=Nmu,use_initial=use_initial)
        rows.append(dict(K=float(K),P_J_S=PS,P_J_H=PH,chiJ=(PH/PS if PS>0 else 0.0)))
        gc.collect()
    return pd.DataFrame(rows)


def initial_parity_control(modefile, K_values, Np=70, Nmu=60):
    return spectrum(modefile,K_values,Np=Np,Nmu=Nmu,use_initial=True)


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--modefile',default='data/raw/linear_handoff_eps1e-3.npz')
    p.add_argument('--out',default='data/derived/relativistic_current_spectrum.csv')
    p.add_argument('--kmin',type=float,default=.2); p.add_argument('--kmax',type=float,default=3.0)
    p.add_argument('--nk',type=int,default=20); p.add_argument('--Np',type=int,default=180); p.add_argument('--Nmu',type=int,default=140)
    return p.parse_args()

if __name__=='__main__':
    a=parse_args(); K=np.linspace(a.kmin,a.kmax,a.nk)
    df=spectrum(a.modefile,K,Np=a.Np,Nmu=a.Nmu)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True); df.to_csv(a.out,index=False)
    ctl=initial_parity_control(a.modefile,K,Np=max(40,a.Np//2),Nmu=max(30,a.Nmu//2))
    ctl.to_csv(Path(a.out).with_name('relativistic_current_parity_control.csv'),index=False)
    print(df.to_string(index=False))
    print('max |initial parity chi|=',ctl.chiJ.abs().max())
