from __future__ import annotations
import json
from dataclasses import asdict
from typing import Any
import numpy as np
from scipy.optimize import minimize_scalar
from fgt_csaudit.pb import SupportLock, scale_xy, _interpolate
from fgt_csaudit.correction_fit import NuFit


def support_from_jsonable(d: dict[str,Any]) -> SupportLock:
    return SupportLock(
      tuple(int(x) for x in d['sizes']),
      {int(k):np.asarray(v,dtype=int) for k,v in d['target_indices'].items()},
      float(d['x_low']),float(d['x_high']),int(d['n_ordered_residuals']),
      float(d['reference_tc']),float(d['reference_nu']),tuple(float(x) for x in d['x_window']),int(d['edge_points'])
    )


def pb_residue_source_masked(curves: dict[int,tuple[np.ndarray,np.ndarray]], *, tc: float, nu: float, channel: str,
                             support: SupportLock, source_ok: dict[int,np.ndarray], interpolation: str='cubic4', invalid_penalty: float=1e6) -> tuple[float,int,bool]:
    if interpolation!='cubic4':
        raise ValueError('strict ESS source-mask sensitivity is locked to the primary cubic4 interpolation')
    if tc<=0 or nu<=0 or not np.isfinite([tc,nu]).all(): return float(invalid_penalty),support.n_ordered_residuals,False
    scaled={}
    for L in support.sizes:
        t,y=curves[L]; x,yy=scale_xy(t,y,L,tc,nu,0.0,channel)
        if np.any(np.diff(x)<=0): return float(invalid_penalty),support.n_ordered_residuals,False
        ok=np.asarray(source_ok[int(L)],dtype=bool)
        if len(ok)!=len(x): raise ValueError(f'L={L}: source mask length mismatch')
        scaled[int(L)]=(x,yy,ok)
    residual_sum=0.0; n=0
    for base_L in support.sizes:
        xb_all,yb_all,source_pass=scaled[int(base_L)]
        base_finite=np.isfinite(xb_all)&np.isfinite(yb_all)&source_pass
        xb=xb_all[base_finite]; yb=yb_all[base_finite]
        if len(xb)<4: return float(invalid_penalty),support.n_ordered_residuals,False
        for target_L in support.sizes:
            if target_L==base_L: continue
            xt_all,yt_all,_=scaled[int(target_L)]
            idx=support.target_indices[int(target_L)]; xt=xt_all[idx]; yt=yt_all[idx]
            if np.any(~np.isfinite(xt)) or np.any(~np.isfinite(yt)): return float(invalid_penalty),support.n_ordered_residuals,False
            if np.any(xt<xb[0]) or np.any(xt>xb[-1]): return float(invalid_penalty),support.n_ordered_residuals,False
            try: pred=_interpolate(xb,yb,xt,'cubic4')
            except Exception: return float(invalid_penalty),support.n_ordered_residuals,False
            if np.any(~np.isfinite(pred)): return float(invalid_penalty),support.n_ordered_residuals,False
            residual_sum += float(np.sum(np.abs(yt-pred))); n += int(len(xt))
    if n != support.n_ordered_residuals: return float(invalid_penalty),n,False
    return residual_sum/float(n),n,True


def fit_nu_source_masked(curves, *, tc: float, channel: str, support: SupportLock, spec: dict[str,Any], nu_bounds: tuple[float,float], source_ok: dict[int,np.ndarray]) -> NuFit:
    lo,hi=map(float,nu_bounds); penalty=float(spec['nu_fit']['invalid_penalty'])
    if not (lo<1.0<hi): raise ValueError('locked symmetric nu bounds must bracket 1')
    grid=np.linspace(lo,hi,int(spec['nu_fit']['coarse_points']))
    def f(nu: float) -> float:
        return float(pb_residue_source_masked(curves,tc=float(tc),nu=float(nu),channel=channel,support=support,source_ok=source_ok,invalid_penalty=penalty)[0])
    vals=np.asarray([f(float(x)) for x in grid],float); valid=np.isfinite(vals)&(vals<0.1*penalty)
    if not np.any(valid): return NuFit(channel,float(tc),float('nan'),penalty,False,True,(lo,hi),len(grid),'no valid coarse points after source mask; no fallback')
    iv=np.flatnonzero(valid); ibest=int(iv[np.argmin(vals[iv])]); best_nu=float(grid[ibest]); best_pb=float(vals[ibest]); nfev=len(grid)
    left=max(0,ibest-1); right=min(len(grid)-1,ibest+1)
    if left<right:
        a,b=float(grid[left]),float(grid[right])
        try:
            res=minimize_scalar(f,bounds=(a,b),method='bounded',options={'xatol':float(spec['nu_fit']['refine_xatol'])})
            nfev += int(getattr(res,'nfev',0)); score=float(f(float(res.x))); nfev+=1
            if bool(res.success) and np.isfinite(score) and score<best_pb: best_nu,best_pb=float(res.x),score
        except Exception: pass
    width=hi-lo; boundary=(best_nu-lo)<=0.02*width or (hi-best_nu)<=0.02*width or ibest in {0,len(grid)-1}
    return NuFit(channel,float(tc),best_nu,best_pb,True,bool(boundary),(lo,hi),int(nfev),'same v3.2.1 optimizer; failed ESS cells excluded only as cubic4 base/source points')


def dependency_cells_at_nu(curves, *, tc: float, nu: float, channel: str, support: SupportLock, source_ok: dict[int,np.ndarray]) -> set[tuple[int,int,float]]:
    """Unique failed source cells selected by the original cubic4 stencil at this (tc,nu).

    Targets remain the locked original targets. This is an attribution diagnostic only.
    """
    scaled={}
    for L in support.sizes:
        t,y=curves[int(L)]; x,yy=scale_xy(t,y,int(L),tc,nu,0.0,channel)
        scaled[int(L)]=(np.asarray(t,float),x,yy,np.asarray(source_ok[int(L)],bool))
    used=set()
    for base_L in support.sizes:
        tbase,xb_all,yb_all,pass_all=scaled[int(base_L)]
        finite=np.isfinite(xb_all)&np.isfinite(yb_all)
        orig_idx=np.flatnonzero(finite); xb=xb_all[finite]
        if len(xb)<4: continue
        for target_L in support.sizes:
            if target_L==base_L: continue
            _,xt_all,yt_all,_=scaled[int(target_L)]; idx=support.target_indices[int(target_L)]
            xt=xt_all[idx]; yt=yt_all[idx]
            if np.any(~np.isfinite(xt)) or np.any(~np.isfinite(yt)): continue
            if np.any(xt<xb[0]) or np.any(xt>xb[-1]): continue
            k=np.searchsorted(xb,xt,side='left'); starts=np.clip(k-2,0,len(xb)-4)
            inds=starts[:,None]+np.arange(4)[None,:]
            for j in np.unique(inds):
                oi=int(orig_idx[int(j)])
                if not bool(pass_all[oi]): used.add((int(base_L),oi,float(tbase[oi])))
    return used
