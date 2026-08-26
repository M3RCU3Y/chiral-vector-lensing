#!/usr/bin/env python3
"""Regenerate the corrected v11 numerical figures from public data.

The public repository tracks generators and numerical inputs rather than PNG
binaries. Analytic manuscript-only figures are outside this calculation repo.
"""
from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
RES=ROOT/'results'/'corrected'; DER=ROOT/'data'/'derived'


def save(fig,name):
    fig.tight_layout()
    fig.savefig(FIG/name,dpi=220,bbox_inches='tight')
    plt.close(fig)


def handoff():
    d=np.load(ROOT/'data/raw/linear_handoff_eps1e-3.npz')
    k=d['kappa']; a=float(d['alpha']); r=float(d['mass_ratio'])
    w=np.sqrt(k*k+(a*r)**2)
    np_=np.maximum(np.asarray(d['occupation_plus'],float),1e-30)
    nm_=np.maximum(np.asarray(d['occupation_minus'],float),1e-30)
    fig,ax=plt.subplots(figsize=(6.6,4.1))
    ax.semilogy(k,np_,label=r'$s=+1$')
    ax.semilogy(k,nm_,label=r'$s=-1$')
    ax.set_xlabel(r'$k/(a_*m_a)$'); ax.set_ylabel('occupation number')
    ax.set_title(rf'Vacuum-to-classical handoff: $\rho_A/\rho_\phi={float(d["epsilon"]):.3e}$')
    ax.legend(frameon=False); ax.grid(alpha=.18)
    save(fig,'linear_to_classical_handoff_corrected.png')


def relativistic_current():
    d=pd.read_csv(DER/'relativistic_current_spectrum.csv')
    ctl=pd.read_csv(DER/'relativistic_current_parity_control.csv')
    fig,ax=plt.subplots(figsize=(6.6,4.1))
    ax.plot(d.K,d.chiJ,marker='o',label='AEN two-helicity handoff')
    ax.plot(ctl.K,ctl.chiJ,linestyle='--',label='parity-symmetric vacuum control')
    ax.axhline(0,linewidth=.8)
    ax.axhline(1,linewidth=.6,linestyle=':'); ax.axhline(-1,linewidth=.6,linestyle=':')
    ax.set_xlabel(r'$K/(a_*m_a)$'); ax.set_ylabel(r'$\chi_J=P_{J,H}/P_{J,S}$')
    ax.set_ylim(-1.05,1.05); ax.grid(alpha=.18); ax.legend(frameon=False)
    ax.set_title('Relativistic Gaussian momentum-current helicity')
    save(fig,'aen_full_two_helicity_current_corrected.png')


def nonlinear_robustness():
    fig,axs=plt.subplots(1,2,figsize=(10.4,3.8),sharey=True)
    for seed in (260823,260824,260825):
        d=pd.read_csv(RES/f'nonlinear_N24_dt0.02_seed{seed}_full.csv')
        axs[0].plot(d.time,d.chiJ_vector,label=str(seed))
    axs[0].axhline(0,linewidth=.7); axs[0].set_title(r'Stochastic seeds, $24^3$')
    axs[0].set_xlabel('simulation time'); axs[0].set_ylabel(r'$\chi_J$'); axs[0].legend(frameon=False,fontsize=8)
    for N in (25,32,40):
        fn=(RES/'nonlinear_N25_dt0.02_seed260823.csv') if N==25 else (RES/f'nonlinear_N{N}_dt0.02_seed260823_full.csv')
        d=pd.read_csv(fn)
        axs[1].plot(d.time,d.chiJ_vector,label=rf'${N}^3$')
    axs[1].axhline(0,linewidth=.7); axs[1].set_title(r'Common-support hierarchy, $k_{\rm init}=2.4$')
    axs[1].set_xlabel('simulation time'); axs[1].legend(frameon=False,fontsize=8)
    for ax in axs: ax.grid(alpha=.18); ax.set_ylim(-.75,.65)
    save(fig,'nonlinear_helicity_seed_robustness_corrected.png')


def backreaction_control():
    nl=pd.read_csv(RES/'nonlinear_N24_dt0.02_seed260823_full.csv')
    li=pd.read_csv(RES/'linear_control_N24_dt0.02_seed260823.csv')
    fig,ax=plt.subplots(figsize=(6.6,4.1))
    ax.plot(nl.time,nl.chiJ_vector,label='nonlinear backreaction')
    ax.plot(li.time,li.chiJ_vector,label='matched no-backreaction control')
    ax.axhline(0,linewidth=.8); ax.set_xlabel('simulation time'); ax.set_ylabel(r'$\chi_J$')
    ax.set_title('Backreaction creates the opposite-sign excursion')
    ax.grid(alpha=.18); ax.legend(frameon=False)
    save(fig,'backreaction_vs_linear_helicity_corrected.png')


def helicity_heatmap():
    d=pd.read_csv(RES/'nonlinear_shells_N24_dt0.01_seed260823_full.csv')
    piv=d.pivot(index='k_mean',columns='time',values='chi_J_shell').sort_index()
    fig,ax=plt.subplots(figsize=(7.0,4.2))
    m=ax.pcolormesh(piv.columns,piv.index,piv.values,shading='auto',vmin=-1,vmax=1,cmap='coolwarm')
    fig.colorbar(m,ax=ax,label=r'$\chi_J(k,t)$')
    ax.set_xlabel('simulation time'); ax.set_ylabel(r'shell $k/(a_*m_a)$')
    ax.set_title('Scale-dependent current-helicity transport')
    save(fig,'helicity_cascade_heatmap_corrected.png')


def kinematic():
    d=pd.read_csv(RES/'nonlinear_N24_dt0.01_seed260823_full.csv')
    fig,ax=plt.subplots(figsize=(5.2,4.1))
    ax.plot(d.time,d.beta_energy_mean,label=r'energy-weighted $\langle p/m_{\gamma\prime}\rangle$')
    ax.axhline(1,linestyle='--',linewidth=.8,label='nonrelativistic threshold')
    ax.set_xlabel('simulation time'); ax.set_ylabel(r'$p/m_{\gamma\prime}$')
    ax.grid(alpha=.18); ax.legend(frameon=False,fontsize=8)
    ax.set_title('The resolved vector population stays relativistic')
    save(fig,'nonlinear_kinematic_evolution_corrected.png')


def iruv():
    d=pd.read_csv(RES/'nonlinear_shells_N24_dt0.01_seed260823_full.csv')
    rows=[]
    for t,g in d.groupby('time'):
        tot=g.vector_mode_energy.sum()
        rows.append((t,g.loc[g.k_mean<1,'vector_mode_energy'].sum()/tot,
                     g.loc[g.k_mean>2.5,'vector_mode_energy'].sum()/tot,
                     (g.k_mean*g.vector_mode_energy).sum()/tot))
    q=pd.DataFrame(rows,columns=['time','ir','uv','kmean'])
    fig,ax=plt.subplots(figsize=(5.2,4.1))
    ax.semilogy(q.time,np.maximum(q.ir,1e-14),label=r'$k<1$')
    ax.semilogy(q.time,np.maximum(q.uv,1e-14),label=r'$k>2.5$')
    ax.set_xlabel('simulation time'); ax.set_ylabel('fraction of shell spectral energy')
    ax.grid(alpha=.18); ax.legend(frameon=False); ax.set_title('Resolved infrared and ultraviolet transfer')
    save(fig,'resolved_IR_UV_transfer_corrected.png')


def dynamic_range():
    r=.5; kpk=1.47; eta=2.
    v=np.logspace(-4,-1,500)
    beta=v/np.sqrt(1-v*v); kcap=r*beta
    R=eta*kpk/kcap
    N=3*np.ceil(R)+1
    fig,ax=plt.subplots(figsize=(6.6,4.1))
    ax.loglog(v,N,label='one Fourier interval below capture')
    for nir,ls in ((4,'--'),):
        ax.loglog(v,3*np.ceil(eta*nir*kpk/kcap)+1,linestyle=ls,label=f'{nir} intervals below capture')
    tab=pd.read_csv(DER/'corrected_capture_dynamic_range.csv')
    for case,g in tab.groupby('case'):
        row=g[g.resolved_IR_modes_below_capture==1].iloc[0]
        ax.scatter([row.escape_speed_over_c],[row.minimum_N_2over3_dealiased],s=35)
        short='pilot' if case.startswith('corrected') else 'AEN displayed halo'
        ax.annotate(short,(row.escape_speed_over_c,row.minimum_N_2over3_dealiased),xytext=(5,5),textcoords='offset points',fontsize=8)
    ax.set_xlabel(r'escape speed $v_{\rm esc}/c$'); ax.set_ylabel('minimum cells per dimension')
    ax.grid(alpha=.18,which='both'); ax.legend(frameon=False,fontsize=8)
    ax.set_title('Strict-dealiased uniform-grid capture requirement')
    save(fig,'capture_dynamic_range_requirement_corrected.png')


def main():
    handoff(); relativistic_current(); nonlinear_robustness(); backreaction_control()
    helicity_heatmap(); kinematic(); iruv(); dynamic_range()
    print('generated corrected manuscript figures in',FIG)

if __name__=='__main__': main()
