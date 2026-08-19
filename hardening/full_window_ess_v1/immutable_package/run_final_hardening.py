from __future__ import annotations
import argparse, json, os, shutil, sys, zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd

PACKAGE_ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(PACKAGE_ROOT/'src'))
sys.path.insert(0,str(PACKAGE_ROOT/'vendor'))

from honeycomb_hardening.common import *
from honeycomb_hardening.ess import build_cell_audit, source_mask_arrays
from honeycomb_hardening.masked_fit import support_from_jsonable, fit_nu_source_masked, dependency_cells_at_nu
from fgt_csaudit.io import build_case
from fgt_csaudit.correction_fit import fit_nu_fixed_tc, restrict_curves

LOCK=load_json(CONFIG_DIR/'INPUT_HASH_LOCK.json')
METHOD=load_json(PACKAGE_ROOT/'METHOD_INTENT_LOCK.json')
SPEC_PATH=PRIMARY_LOCK/'SPEC_LOCK_v321.json'
ASYM_DRAW=ASYM_REF/'tables'/'BOOTSTRAP_DRAW_LEVEL_SYM_VS_ASYM.csv'


def _resolve(cfg: dict[str,Any]) -> tuple[Path,dict[str,Path],Path]:
    root=Path(cfg['project_root']).expanduser().resolve(); paths=project_paths(root)
    out=(PACKAGE_ROOT/cfg.get('output_dir','OUTPUT_FINAL_HARDENING')).resolve()
    return root,paths,out


def _basic_table_shape(path: Path, label: str, expected_rows: int|None=None) -> dict[str,Any]:
    df=pd.read_csv(path)
    required={'case_label','L','temperature','realization','p_target','m2','m4','xi_over_L','ess_energy','tau_energy','measurement_stride'}
    miss=sorted(required-set(df.columns))
    if miss: raise RuntimeError(f'{label}: missing columns {miss}')
    dup=int(df.duplicated(['case_label','L','temperature','realization']).sum())
    if dup: raise RuntimeError(f'{label}: {dup} duplicate scientific rows')
    if expected_rows is not None and len(df)!=expected_rows: raise RuntimeError(f'{label}: expected {expected_rows} rows, got {len(df)}')
    return {'rows':int(len(df)),'columns':int(len(df.columns)),'case_labels':sorted(map(str,df.case_label.unique())),'sizes':sorted(map(int,df.L.unique())),'duplicates':dup}


def stage_precheck(cfg: dict[str,Any], diluted_only: bool=False) -> dict[str,Any]:
    package_integrity=verify_immutable_package_manifest()
    root,paths,out=_resolve(cfg); out.mkdir(parents=True,exist_ok=True)
    exp=LOCK['expected']; checks=[]
    def hcheck(name: str,path: Path,expected: str,required=True):
        exists=path.exists(); got=sha256_file(path) if exists else None; ok=exists and got==expected
        checks.append({'name':name,'path':str(path),'required':required,'exists':exists,'sha256':got,'expected_sha256':expected,'pass':bool(ok if required else (not exists or ok))})
        if required and not ok: raise RuntimeError(f'FAIL-CLOSED: {name} hash/path mismatch: {path}; expected {expected}, got {got}')
    hcheck('fine_realizations_n60.csv',paths['fine'],exp['fine_realizations_n60.csv'])
    require_pristine=bool(cfg.get('strict_require_pristine',True)) and not diluted_only
    hcheck('reference_realizations.csv',paths['reference'],exp['reference_realizations.csv'],required=require_pristine)
    hcheck('LOCKED_SUPPORT_AND_SYMMETRIC_NU_BOUNDS.json',paths['locked_support'],exp['locked_support_json'])
    hcheck('SPEC_LOCK_USED.json',paths['spec_used'],exp['spec_lock_used_json'])
    hcheck('embedded SPEC_LOCK_v321',SPEC_PATH,exp['embedded_SPEC_LOCK_v321'])
    for n,sha in LOCK['asymmetric_completed'].items(): hcheck('asymmetric:'+n,ASYM_REF/n,sha)
    # exact audited sources most directly relevant to this sensitivity
    hcheck('vendored v321 pb.py',PACKAGE_ROOT/'vendor'/'fgt_csaudit'/'pb.py',exp['v321_pb.py'])
    hcheck('vendored v321 correction_fit.py',PACKAGE_ROOT/'vendor'/'fgt_csaudit'/'correction_fit.py',exp['v321_correction_fit.py'])
    hcheck('vendored v321 cli.py',PACKAGE_ROOT/'vendor'/'fgt_csaudit'/'cli.py',exp['v321_cli.py'])
    # semantic/integrity gates
    shape_f=_basic_table_shape(paths['fine'],'fine_realizations_n60',15300)
    if shape_f['case_labels']!=['random_p080','random_p085','random_p090'] or shape_f['sizes']!=[40,60,80,100,120]: raise RuntimeError('FAIL-CLOSED: diluted N60 table topology mismatch')
    shape_r=None
    if paths['reference'].exists():
        shape_r=_basic_table_shape(paths['reference'],'reference_realizations',1050)
        if shape_r['case_labels']!=['pristine_p100'] or shape_r['sizes']!=[40,60,80,100,120]: raise RuntimeError('FAIL-CLOSED: pristine table topology mismatch')
    quality=load_json(paths['quality_gate']) if paths['quality_gate'].exists() else {}
    adaptive=load_json(paths['adaptive']) if paths['adaptive'].exists() else {}
    spec_lock=load_json(SPEC_PATH)
    expected_quality=str(spec_lock['required_project_gates']['quality_gate_status'])
    expected_adaptive=str(spec_lock['required_project_gates']['adaptive_final_status'])
    if str(quality.get('status'))!=expected_quality: raise RuntimeError(f"FAIL-CLOSED: quality_gate status is {quality.get('status')!r}, expected {expected_quality!r}")
    if str(adaptive.get('status'))!=expected_adaptive: raise RuntimeError(f"FAIL-CLOSED: adaptive_final_decision status is {adaptive.get('status')!r}, expected {expected_adaptive!r}")
    replay=pd.read_csv(ASYM_REF/'tables'/'SYMMETRIC_REPLAY_AUDIT.csv')
    if not bool(replay['pass'].all()): raise RuntimeError('FAIL-CLOSED: embedded completed symmetric replay has a failed row')
    draws=pd.read_csv(ASYM_DRAW,usecols=['case_label','bootstrap_index'])
    for lab in ['random_p080','random_p085','random_p090','pristine_p100']:
        g=draws[draws.case_label==lab]
        if len(g)!=2000 or set(g.bootstrap_index.astype(int))!=set(range(2000)): raise RuntimeError(f'FAIL-CLOSED: asymmetric draw table incomplete for {lab}')
    result={'stage':'PRECHECK','status':'PASS','created_utc':utc_now(),'project_root':str(root),'diluted_only_mode':bool(diluted_only),'checks':checks,'fine_shape':shape_f,'reference_shape':shape_r,'quality_gate_status':expected_quality,'adaptive_status':expected_adaptive,'embedded_full_symmetric_replay_rows':int(len(replay)),'package_integrity':package_integrity,'environment':environment_record()}
    atomic_json(out/'manifests'/'PRECHECK_RESULT.json',result); atomic_text(out/'manifests'/'PRECHECK_PASS.txt','PASS\n')
    print('PRECHECK PASS. No original scientific file was modified.')
    return result


def _load_cases(paths: dict[str,Path], labels: list[str]):
    fine=pd.read_csv(paths['fine']); ref=pd.read_csv(paths['reference']) if paths['reference'].exists() else pd.DataFrame()
    out={}
    for lab in labels:
        p=CASE_P[lab]; src=ref if lab=='pristine_p100' else fine
        if src.empty: raise RuntimeError(f'{lab}: required raw table unavailable')
        out[lab]=build_case(src,lab,p,'reference' if lab=='pristine_p100' else 'fine')
    return out


def _support_bundle(paths: dict[str,Path], labels: list[str]):
    d=load_json(paths['locked_support']); bundle={}
    for lab in labels:
        rec=d[lab]; bundle[lab]={'branch_center':rec['branch_center'],'supports':{},'bounds':{}}
        for ch in PRIMARY_CHANNELS:
            bundle[lab]['supports'][ch]=rec['supports']['full'][ch]
            bundle[lab]['bounds'][ch]=rec['nu_bounds']['full'][ch]
    return bundle


def stage_replay(cfg: dict[str,Any], diluted_only: bool=False) -> dict[str,Any]:
    stage_precheck(cfg,diluted_only=diluted_only)
    root,paths,out=_resolve(cfg); labels=[x for x in cfg['cases'] if (x in DILUTED or not diluted_only)]
    cases=_load_cases(paths,labels); spec=load_json(SPEC_PATH); bundle=_support_bundle(paths,labels)
    draws=pd.read_csv(ASYM_DRAW); rows=[]; tol=5e-10
    # Transparent sentinel set: first, midpoint, last. The completed embedded full replay (52 metrics) is verified cryptographically in precheck.
    sentinels=(0,999,1999)
    for lab in labels:
        gd=draws[draws.case_label==lab].set_index('bootstrap_index')
        case=cases[lab]
        for b in sentinels:
            r=gd.loc[b]; curves=case.bootstrap_curves(int(b),int(spec['base_seed'])); tc=float(r.tc_full_joint)
            for ch in PRIMARY_CHANNELS:
                c=restrict_curves(curves[ch],spec['size_windows']['full']); sp=support_from_jsonable(bundle[lab]['supports'][ch]); bounds=tuple(bundle[lab]['bounds'][ch])
                fit=fit_nu_fixed_tc(c,tc=tc,channel=ch,support=sp,spec=spec,nu_bounds=bounds)
                exp_nu=float(r[f'nu_sym_{ch}_full']); exp_pb=float(r[f'pb_sym_{ch}_full'])
                rec={'case_label':lab,'bootstrap_index':b,'channel':ch,'nu_replay':fit.nu,'nu_frozen':exp_nu,'nu_abs_diff':abs(fit.nu-exp_nu),'pb_replay':fit.pb,'pb_frozen':exp_pb,'pb_abs_diff':abs(fit.pb-exp_pb),'boundary_replay':bool(fit.boundary_hit),'boundary_frozen':bool(r[f'boundary_sym_{ch}_full']),'valid_replay':bool(fit.valid),'valid_frozen':bool(r[f'valid_sym_{ch}_full'])}
                rec['pass']=rec['nu_abs_diff']<=tol and rec['pb_abs_diff']<=tol and rec['boundary_replay']==rec['boundary_frozen'] and rec['valid_replay']==rec['valid_frozen']; rows.append(rec)
    df=pd.DataFrame(rows); atomic_csv(out/'tables'/'LOCAL_SENTINEL_SYMMETRIC_REPLAY.csv',df)
    if not bool(df['pass'].all()): raise RuntimeError('FAIL-CLOSED: local sentinel replay does not reproduce frozen symmetric draws')
    result={'stage':'BASELINE_REPLAY_VERIFICATION','status':'PASS','created_utc':utc_now(),'local_sentinels_per_case':list(sentinels),'local_rows':int(len(df)),'tolerance':tol,'embedded_completed_full_replay':'PASS (verified in precheck)','interpretation':'Local deterministic recreation plus cryptographically verified completed 52-metric full symmetric replay. This stage does not change any scientific result.'}
    atomic_json(out/'manifests'/'BASELINE_REPLAY_RESULT.json',result); print('BASELINE REPLAY VERIFICATION PASS.')
    return result


def stage_ess(cfg: dict[str,Any], diluted_only: bool=False) -> dict[str,Any]:
    stage_precheck(cfg,diluted_only=diluted_only)
    root,paths,out=_resolve(cfg); frames=[pd.read_csv(paths['fine'])]
    if paths['reference'].exists() and not diluted_only: frames.append(pd.read_csv(paths['reference']))
    all_df=pd.concat(frames,ignore_index=True)
    cells=build_cell_audit(all_df,100.0,0.90); atomic_csv(out/'tables'/'ESS_CELL_AUDIT_REBUILT.csv',cells)
    summary=(cells.groupby('case_label').agg(total_cells=('primary_source_cell_pass','size'),failed_energy_cells=('primary_source_cell_pass',lambda s:int((~s).sum())),min_corrected_ESS_E=('essE_min','min'),median_of_cell_median_ESS_E=('essE_median','median')).reset_index())
    atomic_csv(out/'tables'/'ESS_CELL_AUDIT_SUMMARY.csv',summary)
    expected=pd.read_csv(PRIMARY_LOCK/'EXPECTED_DILUTED_FULL_WINDOW_ESS_SUMMARY.csv')
    got=summary[summary.case_label.isin(DILUTED)][['case_label','total_cells','failed_energy_cells']].sort_values('case_label').reset_index(drop=True)
    expc=expected[['case_label','total_cells','failed_energy_cells']].sort_values('case_label').reset_index(drop=True)
    expected_counts_match=bool(got.equals(expc))
    if not expected_counts_match: raise RuntimeError(f'FAIL-CLOSED: rebuilt diluted full-window ESS cell counts do not match the frozen audit checkpoint.\nExpected:\n{expc}\nGot:\n{got}')
    result={'stage':'FULL_WINDOW_ESS_AUDIT','status':'PASS','created_utc':utc_now(),'rule_scope':'post-hoc sensitivity diagnostic only; NOT the original near-Tc publication gate','realization_threshold':100.0,'cell_required_pass_fraction':0.90,'primary_source_mask_metric':'corrected energy ESS only','frozen_diluted_count_checkpoint_match':expected_counts_match,'summary':summary.to_dict(orient='records')}
    atomic_json(out/'manifests'/'ESS_AUDIT_RESULT.json',result); print(summary.to_string(index=False)); return result

# multiprocessing globals
_W={}
def _worker_init(project_root: str, labels: list[str], cell_csv: str, support_json: str, spec_json: str):
    global _W
    paths=project_paths(Path(project_root)); cases=_load_cases(paths,labels); cells=pd.read_csv(cell_csv); support=load_json(Path(support_json)); spec=load_json(Path(spec_json))
    _W={'cases':cases,'cells':cells,'support':support,'spec':spec}

def _worker_chunk(task: dict[str,Any]) -> dict[str,Any]:
    lab=task['case_label']; rows=task['draw_rows']; case=_W['cases'][lab]; spec=_W['spec']; cells=_W['cells']; rec=_W['support'][lab]
    source_ok=source_mask_arrays(case,cells); out=[]
    for r in rows:
        b=int(r['bootstrap_index']); curves=case.bootstrap_curves(b,int(spec['base_seed'])); tc=float(r['tc_full_joint'])
        for ch in PRIMARY_CHANNELS:
            c=restrict_curves(curves[ch],spec['size_windows']['full']); sp=support_from_jsonable(rec['supports']['full'][ch]); bounds=tuple(float(x) for x in rec['nu_bounds']['full'][ch])
            baseline_nu=float(r[f'nu_sym_{ch}_full']); used=dependency_cells_at_nu(c,tc=tc,nu=baseline_nu,channel=ch,support=sp,source_ok=source_ok)
            mf=fit_nu_source_masked(c,tc=tc,channel=ch,support=sp,spec=spec,nu_bounds=bounds,source_ok=source_ok)
            out.append({'case_label':lab,'p':float(case.p),'bootstrap_index':b,'channel':ch,'tc_frozen':tc,
              'nu_baseline':baseline_nu,'pb_baseline':float(r[f'pb_sym_{ch}_full']),'boundary_baseline':bool(r[f'boundary_sym_{ch}_full']),'valid_baseline':bool(r[f'valid_sym_{ch}_full']),
              'nu_ess_masked':float(mf.nu),'pb_ess_masked':float(mf.pb),'boundary_ess_masked':bool(mf.boundary_hit),'valid_ess_masked':bool(mf.valid),
              'paired_delta_nu':float(mf.nu-baseline_nu) if mf.valid and np.isfinite(baseline_nu) else float('nan'),
              'dependency_affected':bool(used),'n_failed_source_cells_used':int(len(used)),
              'failed_source_cells_json':json.dumps([{'L':L,'temperature_index':i,'temperature':T} for L,i,T in sorted(used)],separators=(',',':'))})
    return {'task_id':task['task_id'],'rows':out}


def _target_failure_table(cases, cells, support_json) -> pd.DataFrame:
    rows=[]
    for lab,case in cases.items():
        g=cells[cells.case_label==lab]
        lookup={(int(r.L),float(r.temperature)):bool(r.primary_source_cell_pass) for r in g.itertuples()}
        rec=support_json[lab]
        for ch in PRIMARY_CHANNELS:
            sp=support_from_jsonable(rec['supports']['full'][ch])
            for L in sp.sizes:
                temps=case.temperatures[int(L)]
                for idx in sp.target_indices[int(L)]:
                    T=float(temps[int(idx)]); ok=lookup[(int(L),T)]
                    rows.append({'case_label':lab,'channel':ch,'L':int(L),'temperature_index':int(idx),'temperature':T,'energy_cell_pass':ok,'failed_locked_target':not ok})
    return pd.DataFrame(rows)


def stage_sensitivity(cfg: dict[str,Any], diluted_only: bool=False, force: bool=False) -> dict[str,Any]:
    stage_replay(cfg,diluted_only=diluted_only); stage_ess(cfg,diluted_only=diluted_only)
    root,paths,out=_resolve(cfg); labels=[x for x in cfg['cases'] if (x in DILUTED or not diluted_only)]
    draws=pd.read_csv(ASYM_DRAW); draws=draws[draws.case_label.isin(labels)].copy(); cells=pd.read_csv(out/'tables'/'ESS_CELL_AUDIT_REBUILT.csv'); support_json=load_json(paths['locked_support']); spec=load_json(SPEC_PATH)
    cases=_load_cases(paths,labels)
    tf=_target_failure_table(cases,cells,support_json); atomic_csv(out/'tables'/'LOCKED_TARGET_ESS_AUDIT.csv',tf)
    expected_target=pd.read_csv(PRIMARY_LOCK/'EXPECTED_DILUTED_LOCKED_TARGET_ESS_SUMMARY.csv')
    # Validate each primary channel separately against the frozen target-cell count.
    # This prevents accidental double-counting of identical Binder/xi target maps.
    target_summary=(tf[tf.case_label.isin(DILUTED)].groupby(['case_label','channel']).agg(n_locked_target_cells=('failed_locked_target','size'),failed_locked_target_cells=('failed_locked_target','sum')).reset_index())
    target_summary['failed_locked_target_cells']=target_summary['failed_locked_target_cells'].astype(int)
    e=expected_target[['case_label','n_locked_target_cells','failed_locked_target_cells']].copy(); e['failed_locked_target_cells']=e['failed_locked_target_cells'].astype(int)
    e=e.assign(_k=1).merge(pd.DataFrame({'channel':list(PRIMARY_CHANNELS),'_k':[1,1]}),on='_k').drop(columns='_k')
    e=e[['case_label','channel','n_locked_target_cells','failed_locked_target_cells']]
    if not target_summary.sort_values(['case_label','channel']).reset_index(drop=True).equals(e.sort_values(['case_label','channel']).reset_index(drop=True)):
        raise RuntimeError(f'FAIL-CLOSED: locked-target ESS audit no longer matches the frozen diluted checkpoint.\nExpected:\n{e}\nGot:\n{target_summary}')
    # checkpointed chunks; signature-specific directory prevents stale/mixed runs.
    signature=stable_hash({'method':METHOD['ess_dependency_sensitivity'],'fine':sha256_file(paths['fine']),'reference':sha256_file(paths['reference']) if paths['reference'].exists() else None,'locked_support':sha256_file(paths['locked_support']),'draws':sha256_file(ASYM_DRAW),'spec':sha256_file(SPEC_PATH),'labels':labels,'immutable_package_manifest':sha256_file(PACKAGE_ROOT/'IMMUTABLE_PACKAGE_SHA256_MANIFEST.csv'),'runner':sha256_file(PACKAGE_ROOT/'run_final_hardening.py'),'ess_module':sha256_file(PACKAGE_ROOT/'src'/'honeycomb_hardening'/'ess.py'),'masked_fit_module':sha256_file(PACKAGE_ROOT/'src'/'honeycomb_hardening'/'masked_fit.py')})
    cp=out/'checkpoints'/'ess_source_mask'/signature[:16]; cp.mkdir(parents=True,exist_ok=True)
    tasks=[]; chunk=25; task_id=0
    for lab in labels:
        gd=draws[draws.case_label==lab].sort_values('bootstrap_index')
        records=gd.to_dict(orient='records')
        for start in range(0,len(records),chunk):
            part=records[start:start+chunk]; stem=f'{lab}_{start:04d}_{start+len(part)-1:04d}'; csvp=cp/(stem+'.csv'); metap=cp/(stem+'.json')
            if not force and csvp.exists() and metap.exists():
                try:
                    m=load_json(metap)
                    if m.get('signature')==signature and m.get('sha256')==sha256_file(csvp): continue
                except Exception: pass
            tasks.append({'task_id':task_id,'case_label':lab,'draw_rows':part,'csv':str(csvp),'meta':str(metap)}); task_id+=1
    workers=max(1,int(cfg.get('workers',4))); print(f'[ESS source-mask] workers={workers} tasks={len(tasks)}',flush=True)
    if tasks:
        with ProcessPoolExecutor(max_workers=workers,initializer=_worker_init,initargs=(str(root),labels,str(out/'tables'/'ESS_CELL_AUDIT_REBUILT.csv'),str(paths['locked_support']),str(SPEC_PATH))) as ex:
            futures={ex.submit(_worker_chunk,t):t for t in tasks}
            done=0
            for fut in as_completed(futures):
                t=futures[fut]; result=fut.result(); df=pd.DataFrame(result['rows']); csvp=Path(t['csv']); metap=Path(t['meta']); atomic_csv(csvp,df); atomic_json(metap,{'signature':signature,'sha256':sha256_file(csvp),'task_id':result['task_id']}); done+=1
                if done%10==0 or done==len(tasks): print(f'[ESS source-mask] complete {done}/{len(tasks)}',flush=True)
    files=sorted(cp.glob('*.csv')); full=pd.concat([pd.read_csv(p) for p in files],ignore_index=True); full=full.sort_values(['case_label','bootstrap_index','channel']).reset_index(drop=True)
    expected=len(labels)*2000*2
    if len(full)!=expected: raise RuntimeError(f'FAIL-CLOSED: expected {expected} draw-channel rows, collected {len(full)}')
    if full.duplicated(['case_label','bootstrap_index','channel']).any(): raise RuntimeError('FAIL-CLOSED: duplicate draw rows after checkpoints')
    atomic_csv(out/'tables'/'ESS_MASKED_DRAW_LEVEL.csv',full)
    # aggregate dependency map by parsing only unique used cells per draw
    dep_counter=Counter(); dep_draw_counter=Counter()
    for r in full.itertuples():
        cells_used=json.loads(r.failed_source_cells_json)
        for c in cells_used:
            key=(r.case_label,r.channel,int(c['L']),int(c['temperature_index']),float(c['temperature']))
            dep_counter[key]+=1; dep_draw_counter[key]+=1
    dep_rows=[{'case_label':k[0],'channel':k[1],'L':k[2],'temperature_index':k[3],'temperature':k[4],'affected_draws':v,'total_draws':2000,'affected_draw_fraction':v/2000.0} for k,v in sorted(dep_draw_counter.items())]
    dep=pd.DataFrame(dep_rows,columns=['case_label','channel','L','temperature_index','temperature','affected_draws','total_draws','affected_draw_fraction']); atomic_csv(out/'tables'/'INTERPOLATION_DEPENDENCY_MAP.csv',dep)
    # summaries
    rows=[]
    for (lab,ch),g in full.groupby(['case_label','channel'],sort=True):
        base=g.nu_baseline.to_numpy(float); mask=g.nu_ess_masked.to_numpy(float); valid=g.valid_ess_masked.astype(bool).to_numpy(); finite=valid&np.isfinite(mask)
        ci=(float(np.nanpercentile(mask[finite],2.5)),float(np.nanpercentile(mask[finite],97.5))) if finite.any() else (float('nan'),float('nan'))
        base_valid=g.valid_baseline.astype(bool).to_numpy(); base_finite=base_valid & np.isfinite(base)
        base_ci=(float(np.nanpercentile(base[base_finite],2.5)),float(np.nanpercentile(base[base_finite],97.5))) if base_finite.any() else (float('nan'),float('nan'))
        d=g.paired_delta_nu.to_numpy(float); d=d[np.isfinite(d)]; dci=(float(np.percentile(d,2.5)),float(np.percentile(d,97.5))) if len(d) else (float('nan'),float('nan'))
        boundary=float(g.boundary_ess_masked.astype(bool).mean()); valid_frac=float(valid.mean()); contains1=bool(np.isfinite(ci).all() and ci[0]<=1.0<=ci[1]); gate=boundary<=0.10+1e-15 and valid_frac>=0.95-1e-15
        base_boundary=float(g.boundary_baseline.astype(bool).mean()); base_valid_frac=float(base_valid.mean()); base_contains1=bool(np.isfinite(base_ci).all() and base_ci[0]<=1.0<=base_ci[1])
        tfailed=int(tf[(tf.case_label==lab)&(tf.channel==ch)].failed_locked_target.sum())
        rows.append({'case_label':lab,'p':float(g.p.iloc[0]),'channel':ch,'n_bootstrap':len(g),'baseline_nu_median':float(np.median(base[base_finite])) if base_finite.any() else float('nan'),'baseline_ci_low':base_ci[0],'baseline_ci_high':base_ci[1],'baseline_boundary_fraction':base_boundary,'baseline_valid_fraction':base_valid_frac,'baseline_ci_contains_1':base_contains1,'masked_nu_median':float(np.nanmedian(mask[finite])) if finite.any() else float('nan'),'masked_ci_low':ci[0],'masked_ci_high':ci[1],'masked_boundary_fraction':boundary,'masked_valid_fraction':valid_frac,'masked_identifiability_gate_pass':gate,'masked_ci_contains_1':contains1,'paired_delta_median':float(np.median(d)) if len(d) else float('nan'),'paired_delta_ci_low':dci[0],'paired_delta_ci_high':dci[1],'dependency_affected_draw_fraction':float(g.dependency_affected.astype(bool).mean()),'max_failed_source_cells_used_in_one_draw':int(g.n_failed_source_cells_used.max()),'failed_locked_target_cells':tfailed})
    summ=pd.DataFrame(rows); atomic_csv(out/'tables'/'ESS_MASKED_NU_RESULTS.csv',summ)
    # case-level interpretation with existing locked gates only; no invented effect-size threshold.
    cases_out=[]
    for lab,g in summ.groupby('case_label',sort=True):
        both_gate=bool(g.masked_identifiability_gate_pass.all()); lows=g.masked_ci_low.to_numpy(float); highs=g.masked_ci_high.to_numpy(float)
        both_exclude_same_side=bool(both_gate and ((lows>1).all() or (highs<1).all()))
        target_fail=int(g.failed_locked_target_cells.max())
        if both_exclude_same_side: status='POTENTIAL_DECISION_IMPACT_REQUIRES_FULL_REVIEW'
        elif target_fail>0: status='TARGET_ESS_WARNING_CONSERVATIVE_CONCLUSION_NOT_OVERTURNED'
        elif not both_gate: status='SOURCE_MASK_REMAINS_NONIDENTIFIABLE_OR_LIMITED_RANGE'
        else: status='NO_JOINT_EVIDENCE_AGAINST_NU1_UNDER_SOURCE_MASK'
        cases_out.append({'case_label':lab,'status':status,'both_masked_gates_pass':both_gate,'both_masked_intervals_exclude_1_same_side':both_exclude_same_side,'failed_locked_target_cells':target_fail})
    diluted=[x for x in cases_out if x['case_label'] in DILUTED]
    conclusion_stable=not any(x['both_masked_intervals_exclude_1_same_side'] for x in diluted)
    strong_closure=conclusion_stable and all(x['failed_locked_target_cells']==0 for x in diluted) and all(float(summ[(summ.case_label==x['case_label'])].masked_valid_fraction.min())>=0.95-1e-15 for x in diluted)
    result={'stage':'ESS_INFERENCE_DEPENDENCY_SENSITIVITY','status':'PASS','created_utc':utc_now(),'run_signature':signature,'scope':'post-hoc shadow sensitivity; primary symmetric decision unchanged','cases':cases_out,'diluted_primary_conclusion_stable':bool(conclusion_stable),'strong_ess_dependency_closure_supported':bool(strong_closure),'strong_ess_dependency_closure_scope':'diluted cases only','pristine_primary_calibration_redefined':False,'closure_definition':'Strong closure requires: no diluted locked target fails the post-hoc energy-ESS cell rule, masked valid fractions >=0.95, and no diluted case has both primary masked intervals exclude nu=1 on the same side with both locked gates passing.','no_new_monte_carlo':True}
    atomic_json(out/'FINAL_ESS_DEPENDENCY_RESULT.json',result)
    report_lines=[
      '# Final ESS inference-dependency sensitivity report', '',
      '**Scope:** post-hoc shadow sensitivity only. The prespecified symmetric v3.2.1/v3.2.1.1 analysis remains primary.', '',
      '- No new Monte Carlo is run.',
      '- Full-window corrected-energy ESS uses `ESS_corrected = ess_energy * measurement_stride`.',
      '- The 100 / 90% full-window rule is a post-hoc sensitivity diagnostic, not the original locked near-Tc publication gate.',
      '- Failed cells are removed only as cubic4 base/source points. Locked targets, per-draw Tc, bootstrap indices, support, symmetric nu bounds, interpolation type, Pb definition otherwise, and decision gates remain fixed.',
      '- No fallback interpolation is allowed; insufficient/bracketing failures are invalid fits.',
      '- The strong-closure decision is explicitly scoped to the diluted cases. The pristine run is retained only as an ancillary stress diagnostic under this later post-hoc full-window rule and cannot redefine the original exact-limit calibration.', '',
      '## Channel-level results', '',
      summ.to_markdown(index=False), '',
      '## Case-level interpretation', '',
      pd.DataFrame(cases_out).to_markdown(index=False), '',
      f'Diluted primary conclusion stable: **{conclusion_stable}**',
      f'Strong ESS-dependency closure supported: **{strong_closure}**', '',
      'The `strong_ess_dependency_closure_supported` flag is deliberately fail-closed and is defined in `FINAL_ESS_DEPENDENCY_RESULT.json`. It must not be interpreted as proof that the post-hoc full-window ESS rule was part of the original protocol.'
    ]
    atomic_text(out/'FINAL_ESS_DEPENDENCY_REPORT.md','\n'.join(report_lines)+'\n')
    print(json.dumps(result,indent=2)); return result


def stage_manuscript(cfg: dict[str,Any], diluted_only: bool=False) -> dict[str,Any]:
    root,paths,out=_resolve(cfg)
    final=load_json(out/'FINAL_ESS_DEPENDENCY_RESULT.json') if (out/'FINAL_ESS_DEPENDENCY_RESULT.json').exists() else stage_sensitivity(cfg,diluted_only=diluted_only)
    ess=pd.read_csv(out/'tables'/'ESS_MASKED_NU_RESULTS.csv'); asym=pd.read_csv(ASYM_REF/'tables'/'ASYMMETRIC_NU_RESULTS.csv'); shadow=pd.read_csv(ASYM_REF/'tables'/'SHADOW_SENSITIVITY_DECISIONS.csv')
    # S4: asymmetric manuscript table
    s4=asym[['case_label','p','channel','symmetric_nu_median','asymmetric_nu_median','asymmetric_nu_ci_low','asymmetric_nu_ci_high','asymmetric_boundary_fraction','asymmetric_valid_fraction','asymmetric_ci_contains_1']].merge(shadow[['case_label','shadow_status']],on='case_label',how='left')
    atomic_csv(out/'manuscript_support'/'SUPPLEMENTARY_TABLE_S4_ASYMMETRIC.csv',s4)
    s5_all=ess[['case_label','p','channel','baseline_nu_median','baseline_ci_low','baseline_ci_high','baseline_boundary_fraction','baseline_valid_fraction','masked_nu_median','masked_ci_low','masked_ci_high','masked_boundary_fraction','masked_valid_fraction','masked_ci_contains_1','paired_delta_median','paired_delta_ci_low','paired_delta_ci_high','dependency_affected_draw_fraction','failed_locked_target_cells']]
    # The manuscript ESS source-dependency table is deliberately diluted-only. The pristine full-window stress result is retained in the complete audit report, not used to redefine the exact-limit calibration.
    s5=s5_all[s5_all['p']<1.0].copy()
    atomic_csv(out/'manuscript_support'/'SUPPLEMENTARY_TABLE_S5_ESS_DEPENDENCY.csv',s5)
    # One reviewer-facing decision summary across primary, asymmetric, and ESS hardening layers.
    final_cases={x['case_label']:x for x in final['cases']}
    decision_rows=[]
    for lab in sorted(set(asym.case_label)):
        ga=asym[asym.case_label==lab].copy(); pval=float(ga.p.iloc[0]); gates=bool(ga.asymmetric_identifiability_gate_pass.all())
        al=ga.asymmetric_nu_ci_low.to_numpy(float); ah=ga.asymmetric_nu_ci_high.to_numpy(float)
        ajoint=bool(gates and ((al>1).all() or (ah<1).all()))
        decision_rows.append({'case_label':lab,'p':pval,'primary_layer':'SYMMETRIC_V321_V321_1_UNCHANGED','primary_interpretation':'PRISTINE_NU1_COMPATIBLE' if pval==1.0 else 'DILUTED_NU_NOT_IDENTIFIABLE','asymmetric_shadow_status':str(shadow.loc[shadow.case_label==lab,'shadow_status'].iloc[0]),'asymmetric_joint_evidence_against_nu1':ajoint,'ess_source_mask_status':final_cases[lab]['status'],'ess_source_mask_joint_evidence_against_nu1':bool(final_cases[lab]['both_masked_intervals_exclude_1_same_side']),'allowed_final_statement':'PRISTINE_PRIMARY_CALIBRATION_UNCHANGED' if pval==1.0 else 'NO_JOINT_EVIDENCE_AGAINST_NU1; PRIMARY_DECISION_UNCHANGED'})
    atomic_csv(out/'manuscript_support'/'FINAL_HARDENING_DECISION_SUMMARY.csv',pd.DataFrame(decision_rows))
    # Rename/copy completed asymmetric figure for manuscript numbering without altering original evidence.
    for ext in ('pdf','png'):
        src=ASYM_REF/'figures'/f'FIG_S1_ASYMMETRIC_NU_SENSITIVITY.{ext}'
        if src.exists(): shutil.copy2(src,out/'manuscript_support'/f'SUPPLEMENTARY_FIGURE_S3_ASYMMETRIC.{ext}')
    stable=bool(final['diluted_primary_conclusion_stable']); strong=bool(final['strong_ess_dependency_closure_supported']); ready=stable and strong
    if ready:
        ess_sentence='The post-hoc ESS source-dependency sensitivity did not overturn the conservative diluted-exponent conclusion. Failed full-window energy-ESS cells were prohibited only as cubic-interpolation source points, while the locked target support, per-draw Tc values, bootstrap indices, symmetric nu bounds, interpolation rule, and decision thresholds were unchanged. No diluted case produced joint Binder and xi/L evidence against nu=1 under this source mask.'
    elif stable:
        ess_sentence='The post-hoc ESS source-dependency sensitivity did not create joint Binder and xi/L evidence against nu=1, but the package-level closure conditions were not all satisfied. This result must be reported as an unresolved sensitivity rather than as closure of the ESS-dependency question.'
    else:
        ess_sentence='The post-hoc ESS source-dependency sensitivity produced a potential decision impact and therefore requires full scientific review before manuscript integration; the primary analysis has not been overwritten.'
    md=f"""# Manuscript hardening insertion plan\n\nStatus: **{'READY FOR CONTROLLED MANUSCRIPT INTEGRATION' if ready else 'HOLD FOR SCIENTIFIC REVIEW'}**\n\nThis file does not overwrite the manuscript. The prespecified symmetric v3.2.1/v3.2.1.1 analysis remains primary.\n\n## 1. Abstract\nReplace the unqualified diluted-boundary sentence with wording equivalent to:\n\n> Under the prespecified symmetric primary analysis, at every p < 1 at least one primary nu channel exceeds the boundary-hit threshold, and all six diluted upper bootstrap endpoints coincide with the symmetric feasible bound. An additional rule-locked asymmetric-domain sensitivity leaves the central estimates essentially unchanged; although boundary fractions fall below the diagnostic threshold at p = 0.85 and 0.90, Binder and xi/L never jointly exclude nu = 1.\n\nDo **not** describe the asymmetric test as prespecified in the original study. The completed audit is a later rule-locked shadow sensitivity.\n\n## 2. End of Sec. 4.3 — asymmetric method paragraph\n> After completion of the primary symmetric-domain analysis, an additional rule-locked post-processing sensitivity analysis was specified before execution. This analysis removed only the forced symmetry of the centrally locked feasible nu domain about nu=1, while leaving the Monte Carlo data, Tc estimator, bootstrap draws, common-support definition, interpolation rule, collapse residue, lattice-size windows, declared search range 0.55<=nu<=1.45, and decision thresholds unchanged. The resulting outputs were treated as shadow sensitivity diagnostics and were not allowed to overwrite the primary v3.2.1/v3.2.1.1 decision layer.\n\n## 3. Sec. 6.4 — asymmetric result\nUse the generated Supplementary Table S4. The correct count is **5 of 6** diluted asymmetric upper endpoints at 1.45; the exception is p=0.90 Binder with upper endpoint 1.273515. Do not translate REQUIRES_FULL_DECISION_REVIEW into an identified exponent.\n\n## 4. Sec. 6.6 — ESS wording\nThe current near-Tc numbers should be explicitly scoped as **the locked near-Tc probe temperature at L=120**. Do not present them as minima over the full simulated window.\n\nAdd the following only as a post-hoc sensitivity result, not an original gate:\n\n> {ess_sentence}\n\n### Diluted ESS source-mask result to report
Use Supplementary Table S5 for exact values. In the completed run, p=0.80 is numerically unchanged because no failed full-window ESS cell is used by the original cubic4 stencils. At p=0.85, the masked analysis remains limited by the locked boundary criterion. At p=0.90, both masked boundary gates pass, while both masked 95% intervals contain nu=1; therefore the source mask does not create joint evidence against nu=1.

The complete audit also retains a pristine post-hoc stress diagnostic. Under the later full-window mask, one pristine locked target cell fails that post-hoc cell rule and the masked pristine collapse becomes invalid. Do not use this later stress mask to redefine or reject the original pristine calibration: the full-window rule was not the original near-Tc gate. Keep this distinction explicit in repository documentation.

## 5. Sec. 7.2 and Conclusions\nWhere the manuscript says “Under the prespecified criterion”, change this to **“Under the prespecified symmetric primary criterion”**. State that the asymmetric-domain and ESS source-dependency checks are additional sensitivities and do not replace the primary decision layer.\n\n## 6. Supplement numbering\n- S5. Asymmetric-domain sensitivity\n- Supplementary Table S4: generated `SUPPLEMENTARY_TABLE_S4_ASYMMETRIC.csv`\n- Supplementary Figure S3: generated copy of the completed asymmetric figure\n- S6. ESS-inference dependency sensitivity\n- Supplementary Table S5: generated `SUPPLEMENTARY_TABLE_S5_ESS_DEPENDENCY.csv` (diluted cases only)
- Repository/reviewer decision summary: `FINAL_HARDENING_DECISION_SUMMARY.csv`\n\n## 7. Submission-only items still required\nAuthors/affiliations/corresponding author and the public repository DOI or persistent URL remain submission metadata and are not invented by this package.\n"""
    atomic_text(out/'manuscript_support'/'MANUSCRIPT_HARDENING_PLAN.md',md)
    guard={'stage':'MANUSCRIPT_SUPPORT','status':'READY_FOR_CONTROLLED_INTEGRATION' if ready else 'HOLD_FOR_SCIENTIFIC_REVIEW','created_utc':utc_now(),'primary_analysis_overwritten':False,'asymmetric_count_lock':'5/6 diluted asymmetric upper endpoints reach 1.45; p0.90 Binder exception 1.2735154573670118','ess_conclusion_stable':stable,'strong_ess_dependency_closure_supported':strong}
    atomic_json(out/'manuscript_support'/'MANUSCRIPT_INTEGRATION_GUARD.json',guard); print(md); return guard


def stage_pack(cfg: dict[str,Any]) -> Path:
    root,paths,out=_resolve(cfg)
    if not (out/'manuscript_support'/'MANUSCRIPT_INTEGRATION_GUARD.json').exists(): raise RuntimeError('Run manuscript-support stage first')
    # The output manifest intentionally excludes itself; self-hashing would become stale immediately after rewrite.
    payload=sorted(x for x in out.rglob('*') if x.is_file() and 'checkpoints' not in x.parts and x.name!='OUTPUT_SHA256_MANIFEST.csv')
    rows=[{'path':p.relative_to(out).as_posix(),'sha256':sha256_file(p),'bytes':p.stat().st_size} for p in payload]
    manifest=pd.DataFrame(rows); atomic_csv(out/'OUTPUT_SHA256_MANIFEST.csv',manifest)
    zip_path=PACKAGE_ROOT/'HONEYCOMB_FINAL_HARDENING_REVIEW_OUTPUT.zip'
    if zip_path.exists(): zip_path.unlink()
    with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in payload: z.write(p,p.relative_to(out).as_posix())
        z.write(out/'OUTPUT_SHA256_MANIFEST.csv','OUTPUT_SHA256_MANIFEST.csv')
    print(f'Packed review output: {zip_path}')
    return zip_path


def main():
    ap=argparse.ArgumentParser(description='Strict final hardening for the quenched site-diluted honeycomb Ising study. No Monte Carlo.')
    ap.add_argument('stage',choices=['precheck','replay','ess','sensitivity','manuscript','pack','all'])
    ap.add_argument('--config',default=str(PACKAGE_ROOT/'USER_CONFIG.json'))
    ap.add_argument('--diluted-only',action='store_true',help='Developer/self-test mode only. Publication run should omit this flag and require pristine raw reference data.')
    ap.add_argument('--force',action='store_true',help='Rebuild ESS sensitivity checkpoints only; never changes original project files.')
    ns=ap.parse_args(); cfg=load_config(Path(ns.config))
    try:
        if ns.stage=='precheck': stage_precheck(cfg,ns.diluted_only)
        elif ns.stage=='replay': stage_replay(cfg,ns.diluted_only)
        elif ns.stage=='ess': stage_ess(cfg,ns.diluted_only)
        elif ns.stage=='sensitivity': stage_sensitivity(cfg,ns.diluted_only,ns.force)
        elif ns.stage=='manuscript': stage_manuscript(cfg,ns.diluted_only)
        elif ns.stage=='pack': stage_pack(cfg)
        elif ns.stage=='all':
            stage_sensitivity(cfg,ns.diluted_only,ns.force); stage_manuscript(cfg,ns.diluted_only); stage_pack(cfg)
    except Exception as exc:
        print(f'FAIL-CLOSED: {type(exc).__name__}: {exc}',file=sys.stderr); raise

if __name__=='__main__': main()
