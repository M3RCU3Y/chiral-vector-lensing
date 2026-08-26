#!/usr/bin/env python3
"""Convergence table for the relativistic Gaussian current quadrature."""
from pathlib import Path
import sys
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from relativistic_current_spectrum import current_covariance

MODE=ROOT/'data/raw/linear_handoff_eps1e-3.npz'
OUT=ROOT/'data/derived/relativistic_current_quadrature_convergence.csv'


def main():
    rows=[]
    for Np,Nmu in ((280,210),(320,240),(400,300)):
        for K in (0.75,1.5,2.25):
            PS,PH,_=current_covariance(MODE,K,Np=Np,Nmu=Nmu)
            rows.append({'Np':Np,'Nmu':Nmu,'K':K,'P_J_S':PS,'P_J_H':PH,'chiJ':PH/PS})
    df=pd.DataFrame(rows)
    ref=df[df.Np==400].set_index('K').chiJ
    df['abs_delta_vs_400']=df.apply(lambda r: abs(r.chiJ-ref.loc[r.K]),axis=1)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(OUT,index=False)
    print(df.to_string(index=False))
    print('max |delta chi| 320x240 vs 400x300 =',df[df.Np==320].abs_delta_vs_400.max())

if __name__=='__main__': main()
