import json,sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'vendor'))
from honeycomb_hardening.common import sha256_file, load_json, verify_immutable_package_manifest
from honeycomb_hardening.ess import build_cell_audit
from honeycomb_hardening.masked_fit import pb_residue_source_masked
from fgt_csaudit.pb import SupportLock, pb_residue

def test_spec_hash():
    lock=load_json(ROOT/'configs'/'INPUT_HASH_LOCK.json')
    assert sha256_file(ROOT/'reference'/'primary_lock'/'SPEC_LOCK_v321.json')==lock['expected']['embedded_SPEC_LOCK_v321']

def test_asym_5_of_6_lock():
    d=pd.read_csv(ROOT/'reference'/'asymmetric_completed'/'tables'/'ASYMMETRIC_NU_RESULTS.csv')
    g=d[d.p<1]
    assert len(g)==6
    assert int(np.isclose(g.asymmetric_nu_ci_high,1.45,atol=1e-12).sum())==5
    ex=g[~np.isclose(g.asymmetric_nu_ci_high,1.45,atol=1e-12)]
    assert len(ex)==1 and ex.iloc[0].case_label=='random_p090' and ex.iloc[0].channel=='binder_roa'

def test_ess_stride_correction_and_rule():
    df=pd.DataFrame({'case_label':['x']*2,'L':[40]*2,'temperature':[1.0]*2,'realization':[0,1],'measurement_stride':[2,2],'ess_energy':[50,49],'ess_abs_m':[100,100],'tau_energy':[1,1],'tau_abs_m':[1,1]})
    a=build_cell_audit(df,100,.9).iloc[0]
    assert a.essE_min==98
    assert a.energy_pass_fraction==0.5
    assert not bool(a.primary_source_cell_pass)

def test_identity_mask_equals_original_pb():
    t=np.linspace(.8,1.2,9)
    curves={40:(t,np.sin(t)),60:(t,np.sin(t)+.01),80:(t,np.sin(t)-.01)}
    target={40:np.array([3,4,5]),60:np.array([3,4,5]),80:np.array([3,4,5])}
    sp=SupportLock((40,60,80),target,-1,1,(3-1)*9,1.0,1.0,(-3,3),2)
    ok={L:np.ones(len(t),bool) for L in curves}
    a=pb_residue(curves,tc=1,nu=1,q=0,channel='binder_roa',support=sp,interpolation='cubic4')
    b=pb_residue_source_masked(curves,tc=1,nu=1,channel='binder_roa',support=sp,source_ok=ok)
    assert a[2] and b[2] and abs(a[0]-b[0])<1e-15

def test_no_fallback_when_mask_leaves_lt4_sources():
    t=np.linspace(.8,1.2,9); curves={40:(t,t),60:(t,t),80:(t,t)}
    target={40:np.array([3,4,5]),60:np.array([3,4,5]),80:np.array([3,4,5])}
    sp=SupportLock((40,60,80),target,-1,1,18,1,1,(-3,3),2)
    ok={L:np.array([1,1,1,0,0,0,0,0,0],bool) for L in curves}
    val,n,valid=pb_residue_source_masked(curves,tc=1,nu=1,channel='binder_roa',support=sp,source_ok=ok)
    assert not valid and val==1e6

def test_frozen_diluted_ess_checkpoints():
    a=pd.read_csv(ROOT/'reference'/'primary_lock'/'EXPECTED_DILUTED_FULL_WINDOW_ESS_SUMMARY.csv').set_index('case_label')
    assert a.loc['random_p080','failed_energy_cells']==13
    assert a.loc['random_p085','failed_energy_cells']==13
    assert a.loc['random_p090','failed_energy_cells']==20
    t=pd.read_csv(ROOT/'reference'/'primary_lock'/'EXPECTED_DILUTED_LOCKED_TARGET_ESS_SUMMARY.csv').set_index('case_label')
    assert t['failed_locked_target_cells'].sum()==0
    assert t.loc['random_p080','n_locked_target_cells']==30
    assert t.loc['random_p085','n_locked_target_cells']==36
    assert t.loc['random_p090','n_locked_target_cells']==36


def test_immutable_package_manifest():
    r=verify_immutable_package_manifest()
    assert r['status']=='PASS'
    assert r['files_verified']>20

def test_source_masked_pb_against_independent_lagrange4():
    # Independent manual four-point Lagrange interpolation for a synthetic masked case.
    t=np.linspace(0.8,1.2,11)
    curves={40:(t,np.cos(1.7*t)),60:(t,np.cos(1.7*t)+0.015),80:(t,np.cos(1.7*t)-0.012)}
    target={40:np.array([3,4,5,6,7]),60:np.array([3,4,5,6,7]),80:np.array([3,4,5,6,7])}
    sp=SupportLock((40,60,80),target,-20,20,30,1.0,1.0,(-3,3),2)
    ok={40:np.ones(11,bool),60:np.ones(11,bool),80:np.ones(11,bool)}
    ok[40][9]=False; ok[60][1]=False; ok[80][9]=False
    got,n,valid=pb_residue_source_masked(curves,tc=1.0,nu=1.0,channel='binder_roa',support=sp,source_ok=ok)
    assert valid and n==30
    def interp4(x,y,xt):
        k=np.searchsorted(x,xt,side='left'); st=np.clip(k-2,0,len(x)-4); out=[]
        for xx,s0 in zip(xt,st):
            xs=x[s0:s0+4]; ys=y[s0:s0+4]; yy=0.0
            for j in range(4):
                w=1.0
                for m in range(4):
                    if m!=j: w*= (xx-xs[m])/(xs[j]-xs[m])
                yy += ys[j]*w
            out.append(yy)
        return np.asarray(out)
    scaled={}
    for L,(tt,yy) in curves.items():
        x=((tt-1.0)/1.0)*float(L)
        scaled[L]=(x,yy)
    total=0.0; nn=0
    for bL in sp.sizes:
        xb,yb=scaled[bL]; m=np.isfinite(xb)&np.isfinite(yb)&ok[bL]; xb=xb[m]; yb=yb[m]
        for tL in sp.sizes:
            if tL==bL: continue
            xt,yt=scaled[tL]; idx=sp.target_indices[tL]; xx=xt[idx]; yy=yt[idx]
            pred=interp4(xb,yb,xx); total += float(np.abs(yy-pred).sum()); nn += len(xx)
    expected=total/nn
    assert abs(got-expected)<1e-13
