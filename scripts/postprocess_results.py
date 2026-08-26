#!/usr/bin/env python3
"""Build the small derived tables used by the manuscript from corrected runs.

This script intentionally performs no expensive field evolution.  It reduces the
checked-in CSV histories into compact convergence/robustness tables, so every
headline number in the paper can be traced to a raw run table.
"""
from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "corrected"
OUT = ROOT / "data" / "derived"
OUT.mkdir(parents=True, exist_ok=True)


def zero_crossings(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    rows = []
    for i in range(len(y)-1):
        if y[i] == 0:
            direction = "neg_to_pos" if y[i+1] > 0 else "pos_to_neg"
            rows.append((float(t[i]), direction))
        elif y[i]*y[i+1] < 0:
            tc = t[i] - y[i]*(t[i+1]-t[i])/(y[i+1]-y[i])
            direction = "neg_to_pos" if y[i] < 0 else "pos_to_neg"
            rows.append((float(tc), direction))
    return rows


def dominant_crossings(df, col="chiJ_vector"):
    t = df["time"].to_numpy(float); y = df[col].to_numpy(float)
    imax = int(np.nanargmax(y)); tpeak = float(t[imax])
    z = zero_crossings(t, y)
    before = [x for x,d in z if x <= tpeak and d == "neg_to_pos"]
    after = [x for x,d in z if x >= tpeak and d == "pos_to_neg"]
    return {
        "first_main_crossing": before[-1] if before else np.nan,
        "second_main_crossing": after[0] if after else np.nan,
        "all_crossings": z,
        "max_positive_chi": float(y[imax]),
        "time_of_max": tpeak,
        "endpoint_chi": float(y[-1]),
    }


def global_energy_drift(df):
    return df["rho_total"].to_numpy(float)/float(df["rho_total"].iloc[0]) - 1.0


def seed_table():
    rows=[]
    for seed in (260823,260824,260825):
        f=RES/f"nonlinear_N24_dt0.02_seed{seed}_full.csv"
        d=pd.read_csv(f)
        a=dominant_crossings(d,"chiJ_vector")
        c=dominant_crossings(d,"chiJ_klt2p5")
        rows.append({
            "seed":seed,
            "first_main_crossing":a["first_main_crossing"],
            "second_main_crossing":a["second_main_crossing"],
            "max_positive_chi":a["max_positive_chi"],
            "time_of_max":a["time_of_max"],
            "endpoint_chi":a["endpoint_chi"],
            "common_k_first":c["first_main_crossing"],
            "common_k_second":c["second_main_crossing"],
            "common_k_peak":c["max_positive_chi"],
            "initial_scale":float(d.initial_scale.iloc[0]),
            "precondition_epsA":float(d.precondition_epsA.iloc[0]),
            "max_lorenz_shadow":float(d.lorenz_shadow_relative_vector_rms.max()),
            "max_abs_energy_drift":float(np.max(np.abs(global_energy_drift(d)))),
        })
    out=pd.DataFrame(rows)
    out.to_csv(OUT/"corrected_seed_robustness.csv",index=False)
    return out


def resolution_table():
    rows=[]
    for N,dt in ((24,.02),(32,.02),(40,.02)):
        f=RES/f"nonlinear_N{N}_dt0.02_seed260823_full.csv"
        d=pd.read_csv(f)
        a=dominant_crossings(d,"chiJ_vector")
        c=dominant_crossings(d,"chiJ_klt2p5")
        rows.append({
            "N":N,"dt":dt,
            "first_main_crossing":a["first_main_crossing"],
            "second_main_crossing":a["second_main_crossing"],
            "max_positive_chi":a["max_positive_chi"],
            "time_of_max":a["time_of_max"],
            "endpoint_chi":a["endpoint_chi"],
            "common_k_first":c["first_main_crossing"],
            "common_k_second":c["second_main_crossing"],
            "common_k_peak":c["max_positive_chi"],
            "endpoint_chi_common":float(d.chiJ_klt2p5.iloc[-1]),
            "n_zero_crossings":len(a["all_crossings"]),
            "max_lorenz_shadow":float(d.lorenz_shadow_relative_vector_rms.max()),
            "max_abs_energy_drift":float(np.max(np.abs(global_energy_drift(d)))),
            "max_longitudinal_fraction":float(d.longitudinal_field_fraction.max()),
            "spherical_safe_kmax":float(d.spherical_safe_kmax.iloc[0]),
            "kinit_effective":float(d.kinit_effective.iloc[0]),
        })
    out=pd.DataFrame(rows)
    out.to_csv(OUT/"corrected_spatial_convergence.csv",index=False)
    return out



def common_support_table():
    """Bridge hierarchy with identical requested kinit=2.4 at N=25,32,40.

    N=25 is the smallest strict-dealiased grid whose isotropic safe sphere
    contains k=2.4.  Its checked-in run extends through t=6.5, which is already
    beyond the second dominant crossing; therefore the table focuses on
    crossing times and peak helicity rather than a common t=7.5 endpoint.
    """
    rows=[]
    files={25:RES/'nonlinear_N25_dt0.02_seed260823.csv',
           32:RES/'nonlinear_N32_dt0.02_seed260823_full.csv',
           40:RES/'nonlinear_N40_dt0.02_seed260823_full.csv'}
    for N,f in files.items():
        d=pd.read_csv(f); a=dominant_crossings(d); c=dominant_crossings(d,'chiJ_klt2p5')
        rows.append({
            'N':N,'dt':float(d.dt.iloc[0]),'t_end':float(d.time.iloc[-1]),
            'kinit_effective':float(d.kinit_effective.iloc[0]),
            'first_main_crossing':a['first_main_crossing'],
            'second_main_crossing':a['second_main_crossing'],
            'max_positive_chi':a['max_positive_chi'],
            'common_k_first':c['first_main_crossing'],
            'common_k_second':c['second_main_crossing'],
            'common_k_peak':c['max_positive_chi'],
            'max_lorenz_shadow':float(d.lorenz_shadow_relative_vector_rms.max()),
        })
    out=pd.DataFrame(rows); out.to_csv(OUT/'corrected_common_support_check.csv',index=False)
    return out

def timestep_table():
    fine=pd.read_csv(RES/"nonlinear_N24_dt0.01_seed260823_full.csv")
    coarse=pd.read_csv(RES/"nonlinear_N24_dt0.02_seed260823_full.csv")
    af=dominant_crossings(fine); ac=dominant_crossings(coarse)
    yc=np.interp(fine.time,coarse.time,coarse.chiJ_vector)
    rows=[{
        "N":24,
        "dt_fine":0.01,"dt_coarse":0.02,
        "fine_first":af["first_main_crossing"],"coarse_first":ac["first_main_crossing"],
        "abs_delta_first":abs(af["first_main_crossing"]-ac["first_main_crossing"]),
        "fine_second":af["second_main_crossing"],"coarse_second":ac["second_main_crossing"],
        "abs_delta_second":abs(af["second_main_crossing"]-ac["second_main_crossing"]),
        "fine_peak":af["max_positive_chi"],"coarse_peak":ac["max_positive_chi"],
        "abs_delta_peak":abs(af["max_positive_chi"]-ac["max_positive_chi"]),
        "fine_endpoint":af["endpoint_chi"],"coarse_endpoint":ac["endpoint_chi"],
        "abs_delta_endpoint":abs(af["endpoint_chi"]-ac["endpoint_chi"]),
        "max_abs_history_delta":float(np.max(np.abs(fine.chiJ_vector-yc))),
        "fine_max_abs_energy_drift":float(np.max(np.abs(global_energy_drift(fine)))),
        "coarse_max_abs_energy_drift":float(np.max(np.abs(global_energy_drift(coarse)))),
    }]
    out=pd.DataFrame(rows)
    out.to_csv(OUT/"corrected_timestep_convergence.csv",index=False)
    return out


def wigner_table():
    conditioned=pd.read_csv(RES/"nonlinear_N24_dt0.02_seed260823_full.csv")
    uncond=pd.read_csv(RES/"nonlinear_N24_dt0.02_seed260823_unconditioned_matched.csv")
    if int(uncond.seed.iloc[0]) != int(conditioned.seed.iloc[0]):
        raise ValueError("Wigner conditioning control must use the same random seed")
    if not np.isclose(float(uncond.initial_scale.iloc[0]), 1.0):
        raise ValueError("unconditioned Wigner control unexpectedly rescaled")
    rows=[]
    for label,d in (("conditioned",conditioned),("unconditioned",uncond)):
        a=dominant_crossings(d)
        rows.append({
            "ensemble":label,
            "initial_scale":float(d.initial_scale.iloc[0]),
            "precondition_epsA":float(d.precondition_epsA.iloc[0]),
            "first_main_crossing":a["first_main_crossing"],
            "second_main_crossing":a["second_main_crossing"],
            "max_positive_chi":a["max_positive_chi"],
            "endpoint_chi":a["endpoint_chi"],
            "max_lorenz_shadow":float(d.lorenz_shadow_relative_vector_rms.max()),
            "max_abs_energy_drift":float(np.max(np.abs(global_energy_drift(d)))),
        })
    out=pd.DataFrame(rows)
    out.to_csv(OUT/"corrected_wigner_conditioning.csv",index=False)
    return out


def capture_table():
    d=pd.read_csv(RES/"nonlinear_N24_dt0.01_seed260823_full.csv")
    vpilot=float(d.escape_speed_over_c.max())
    r=.5; kpk=1.47; eta=2.0
    c_kms=299792.458
    cases=[("corrected pilot endpoint",vpilot),("AEN displayed large-normalization escape speed",537.0/c_kms)]
    rows=[]
    for name,v in cases:
        beta=v/math.sqrt(max(1e-30,1-v*v))
        kcap=r*beta
        for nir in (1,4):
            R=eta*nir*kpk/kcap
            nreq=3*math.ceil(R)+1
            rows.append({
                "case":name,"escape_speed_over_c":v,"k_capture_over_astar_ma":kcap,
                "k_peak":kpk,"UV_margin_times_kpeak":eta,
                "resolved_IR_modes_below_capture":nir,"minimum_N_2over3_dealiased":nreq,
            })
    out=pd.DataFrame(rows)
    out.to_csv(OUT/"corrected_capture_dynamic_range.csv",index=False)
    return out


def validation_table():
    ref=pd.read_csv(RES/"nonlinear_N24_dt0.01_seed260823_full.csv")
    lin=pd.read_csv(RES/"linear_control_N24_dt0.02_seed260823.csv")
    shell=pd.read_csv(RES/"nonlinear_shells_N24_dt0.01_seed260823_full.csv")
    def shell_stats(t):
        s=shell[np.isclose(shell.time,t)]
        tot=float(s.vector_mode_energy.sum())
        return {
            "ir":float(s.loc[s.k_mean<1,"vector_mode_energy"].sum()/tot),
            "uv":float(s.loc[s.k_mean>2.5,"vector_mode_energy"].sum()/tot),
            "kmean":float((s.k_mean*s.vector_mode_energy).sum()/tot),
        }
    s0,sf=shell_stats(0.0),shell_stats(7.5)
    rows=[
        ("reference_max_gauss_relative",float(ref.gauss_constraint_relative_rms.max()),"dimensionless"),
        ("reference_max_lorenz_shadow_relative",float(ref.lorenz_shadow_relative_vector_rms.max()),"dimensionless"),
        ("reference_max_abs_total_energy_drift",float(np.max(np.abs(global_energy_drift(ref)))),"dimensionless"),
        ("reference_max_longitudinal_fraction",float(ref.longitudinal_field_fraction.max()),"dimensionless"),
        ("linear_control_max_longitudinal_fraction",float(lin.longitudinal_field_fraction.max()),"dimensionless"),
        ("linear_control_max_lorenz_shadow_relative",float(lin.lorenz_shadow_relative_vector_rms.max()),"dimensionless"),
        ("reference_beta_energy_initial",float(ref.beta_energy_mean.iloc[0]),"dimensionless"),
        ("reference_beta_energy_final",float(ref.beta_energy_mean.iloc[-1]),"dimensionless"),
        ("reference_max_fraction_beta_lt1",float(ref.frac_beta_lt1.max()),"dimensionless"),
        ("reference_IR_energy_fraction_initial",s0["ir"],"dimensionless"),
        ("reference_IR_energy_fraction_final",sf["ir"],"dimensionless"),
        ("reference_UV_energy_fraction_initial",s0["uv"],"dimensionless"),
        ("reference_UV_energy_fraction_final",sf["uv"],"dimensionless"),
        ("reference_shell_kmean_initial",s0["kmean"],"a_star m_a"),
        ("reference_shell_kmean_final",sf["kmean"],"a_star m_a"),
        ("reference_max_escape_speed_over_c",float(ref.escape_speed_over_c.max()),"dimensionless"),
    ]
    out=pd.DataFrame(rows,columns=["quantity","value","units"])
    out.to_csv(OUT/"corrected_validation_summary.csv",index=False)
    return out


def main():
    for name,fn in [
        ("seed",seed_table),("resolution",resolution_table),("common-support",common_support_table),("timestep",timestep_table),
        ("wigner",wigner_table),("capture",capture_table),("validation",validation_table)
    ]:
        df=fn(); print(f"[{name}] {len(df)} rows -> data/derived")
        print(df.to_string(index=False))

if __name__ == "__main__":
    main()
