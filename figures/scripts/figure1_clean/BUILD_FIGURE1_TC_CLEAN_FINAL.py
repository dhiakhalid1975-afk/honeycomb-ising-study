from pathlib import Path
import hashlib, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'TABLE_TC_RG_CROSSING_BOOTSTRAP.snapshot.csv'
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
EXACT_TC = 1.518651435000414

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

df = pd.read_csv(DATA).sort_values('p')
p = df['p'].to_numpy(float)
med = df['tc_bootstrap_median'].to_numpy(float)
lo = df['tc_bootstrap_ci_low'].to_numpy(float)
hi = df['tc_bootstrap_ci_high'].to_numpy(float)
rl = df['tc_robustness_low'].to_numpy(float)
rh = df['tc_robustness_high'].to_numpy(float)

plt.rcParams.update({
    'font.size': 9.5,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'axes.linewidth': 0.8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
fig, ax = plt.subplots(figsize=(5.9, 4.15))
# Robustness envelope: thick translucent vertical span, explicitly non-probabilistic.
for x,l,h in zip(p,rl,rh):
    ax.plot([x,x],[l,h], lw=9.0, alpha=0.16, solid_capstyle='round', color='0.35', zorder=1)
# 95% bootstrap intervals + median.
yerr=np.vstack([med-lo, hi-med])
ax.errorbar(p, med, yerr=yerr, fmt='o', ms=5.2, mfc='white', mew=1.4,
            capsize=3.0, lw=1.2, color='#24588a', zorder=3)
# Exact pristine point only at p=1.
ax.scatter([1.0],[EXACT_TC], marker='*', s=150, color='#2a6f4e', zorder=4)
ax.set_xlabel(r'Site occupation probability, $p$')
ax.set_ylabel(r'Critical temperature, $T_c/(J/k_B)$')
ax.set_xlim(0.775,1.025)
ax.set_xticks([0.80,0.85,0.90,0.95,1.00])
ax.grid(alpha=0.20)
handles=[
    Line2D([],[],marker='o',mfc='white',mec='#24588a',mew=1.4,ls='none',label=r'joint RG-invariant $T_c$; 95% bootstrap interval'),
    Line2D([],[],color='0.35',lw=8,alpha=0.16,label='robustness envelope (not a confidence interval)'),
    Line2D([],[],marker='*',color='#2a6f4e',ls='none',markersize=10,label='exact pristine honeycomb value'),
]
ax.legend(handles=handles, loc='upper left', framealpha=0.95)
fig.tight_layout()
outputs=[]
for ext in ('pdf','png','svg'):
    path=OUT/f'Figure1_Tc_vs_p_clean.{ext}'
    kwargs={'bbox_inches':'tight'}
    if ext=='png': kwargs['dpi']=600
    fig.savefig(path, **kwargs)
    outputs.append({'path':str(path.name),'sha256':sha256(path),'bytes':path.stat().st_size})
plt.close(fig)
meta={
    'figure':'Figure1_Tc_vs_p_clean',
    'input':'TABLE_TC_RG_CROSSING_BOOTSTRAP.snapshot.csv',
    'input_sha256':sha256(DATA),
    'exact_tc':EXACT_TC,
    'meaning':'Primary dilution-dependent critical temperature. Median and 95% realization-level percentile-bootstrap intervals are shown; the separate robustness envelope is a sensitivity range, not a confidence interval. The exact pristine honeycomb value is an external validation point only. No legacy estimator is plotted.',
    'outputs':outputs
}
with open(OUT/'Figure1_Tc_vs_p_clean.provenance.json','w',encoding='utf-8') as f:
    json.dump(meta,f,indent=2,ensure_ascii=False)
print(json.dumps(meta,indent=2))
