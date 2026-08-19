from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

REQUIRED = {'case_label','L','temperature','realization','measurement_stride','ess_energy','ess_abs_m','tau_energy','tau_abs_m'}


def build_cell_audit(df: pd.DataFrame, realization_threshold: float=100.0, required_pass_fraction: float=0.90) -> pd.DataFrame:
    missing=sorted(REQUIRED-set(df.columns))
    if missing: raise ValueError(f'missing ESS columns: {missing}')
    x=df.copy()
    x['essE_corrected']=pd.to_numeric(x['ess_energy'],errors='coerce')*pd.to_numeric(x['measurement_stride'],errors='coerce')
    x['essM_corrected']=pd.to_numeric(x['ess_abs_m'],errors='coerce')*pd.to_numeric(x['measurement_stride'],errors='coerce')
    x['energy_realization_pass']=x['essE_corrected'] >= float(realization_threshold)
    x['mag_realization_pass']=x['essM_corrected'] >= float(realization_threshold)
    x['joint_realization_pass']=x['energy_realization_pass'] & x['mag_realization_pass']
    def p10(s): return float(np.nanpercentile(np.asarray(s,float),10))
    g=x.groupby(['case_label','L','temperature'],sort=True)
    out=g.agg(
      n=('realization','count'),
      measurement_stride_min=('measurement_stride','min'),measurement_stride_max=('measurement_stride','max'),
      energy_pass_n=('energy_realization_pass','sum'),mag_pass_n=('mag_realization_pass','sum'),joint_pass_n=('joint_realization_pass','sum'),
      essE_min=('essE_corrected','min'),essE_median=('essE_corrected','median'),
      essM_min=('essM_corrected','min'),essM_median=('essM_corrected','median'),
      tauE_max=('tau_energy','max'),tauE_median=('tau_energy','median'),tauM_max=('tau_abs_m','max'),tauM_median=('tau_abs_m','median'),
    ).reset_index()
    p10e=g['essE_corrected'].apply(p10).rename('essE_p10').reset_index()
    p10m=g['essM_corrected'].apply(p10).rename('essM_p10').reset_index()
    out=out.merge(p10e,on=['case_label','L','temperature']).merge(p10m,on=['case_label','L','temperature'])
    out['energy_pass_fraction']=out['energy_pass_n']/out['n']
    out['mag_pass_fraction']=out['mag_pass_n']/out['n']
    out['joint_pass_fraction']=out['joint_pass_n']/out['n']
    out['energy_cell_pass']=out['energy_pass_fraction'] >= float(required_pass_fraction)-1e-15
    out['joint_cell_pass']=out['joint_pass_fraction'] >= float(required_pass_fraction)-1e-15
    out['realization_threshold']=float(realization_threshold)
    out['required_pass_fraction']=float(required_pass_fraction)
    out['primary_mask_metric']='corrected_energy_ESS_only'
    out['primary_source_cell_pass']=out['energy_cell_pass']
    cols=['case_label','L','temperature','n','energy_pass_n','energy_pass_fraction','mag_pass_n','mag_pass_fraction','joint_pass_n','joint_pass_fraction',
          'essE_min','essE_p10','essE_median','essM_min','essM_p10','essM_median','tauE_max','tauE_median','tauM_max','tauM_median',
          'measurement_stride_min','measurement_stride_max','realization_threshold','required_pass_fraction','energy_cell_pass','joint_cell_pass','primary_mask_metric','primary_source_cell_pass']
    return out[cols].sort_values(['case_label','L','temperature']).reset_index(drop=True)


def source_mask_arrays(case, cell_audit: pd.DataFrame) -> dict[int,np.ndarray]:
    g=cell_audit[cell_audit['case_label']==case.label]
    if g.empty: raise ValueError(f'no cell audit rows for {case.label}')
    out={}
    for L in case.sizes:
        h=g[g.L==L].sort_values('temperature')
        lookup={float(r.temperature):bool(r.primary_source_cell_pass) for r in h.itertuples()}
        vals=[]
        for t in case.temperatures[L]:
            ft=float(t)
            if ft not in lookup: raise ValueError(f'{case.label}/L={L}/T={ft}: missing cell audit')
            vals.append(lookup[ft])
        out[int(L)]=np.asarray(vals,dtype=bool)
    return out
