#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BUILD_FIGURE2_rg_crossings.py  (v321_1-fig2-2026-08-12-pairset-comparability)
=========================================================================
Figure 2: representative disorder-averaged RG-invariant curves for the
calibration case (p = 1) and a diluted case, in case-specific plotted windows
centred on each locked Tc, under the same analysis pipeline.

The plotted quantities are the central (disorder-averaged) curves returned by
CaseData.central_curves(): the Binder ratio is formed as a ratio of averages over
quenched realizations, and xi/L is a finite-value mean over realizations. They are
therefore NOT raw Monte Carlo trajectories, and the figure carries no uncertainty
band. Statistical and systematic uncertainty on Tc and nu is reported in the
bootstrap tables and in Figures 3-5, not here.

SCIENTIFIC SCOPE
----------------
This script is a renderer only. It does not estimate, re-estimate, or refine any
critical temperature. The vertical reference line in each panel is the locked joint
RG-crossing Tc read verbatim from TABLE_TC_RG_CROSSING_BOOTSTRAP.csv (column:
tc_bootstrap_median). The legacy fixed-nu susceptibility-shift Tc is never read,
never used as a branch locator, and never enters this figure.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

FIGURE_BUILD_VERSION = "v321_1-fig2-2026-08-12-pairset-comparability"
STEM = "Figure2_rg_invariant_crossings"
FORMATS = ("pdf", "png", "svg")
DPI = 600
ROWS = (("pristine_p100", r"$p=1.00$ (calibration)", 0.032),
        ("random_p080", r"$p=0.80$", 0.030))
COLS = (("binder_roa", r"Binder ratio $U$"), ("xi_over_L", r"$\xi/L$"))

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 8.5, "legend.fontsize": 6.8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.8,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _find(root: Path, name: str) -> Path:
    for c in (root / name, root / "tables" / name):
        if c.exists():
            return c
    hits = sorted(root.rglob(name))
    if not hits:
        raise SystemExit(f"FAIL-CLOSED: required input not found under {root}: {name}")
    return hits[0]


def _adjacent_crossings(curves, sizes, lo, hi, tc, n=4001):
    """Locate adjacent-size crossings inside [lo, hi]. Diagnostic output only."""
    out = []
    for a, b in zip(sizes[:-1], sizes[1:]):
        ta, ya = curves[a]
        tb, yb = curves[b]
        ma = np.isfinite(ta) & np.isfinite(ya)
        mb = np.isfinite(tb) & np.isfinite(yb)
        if ma.sum() < 3 or mb.sum() < 3:
            out.append((int(a), int(b), float("nan"), 0, ""))
            continue
        fa = PchipInterpolator(ta[ma], ya[ma], extrapolate=False)
        fb = PchipInterpolator(tb[mb], yb[mb], extrapolate=False)
        g = np.linspace(max(lo, ta[ma].min(), tb[mb].min()),
                        min(hi, ta[ma].max(), tb[mb].max()), n)
        d = np.asarray(fa(g) - fb(g), float)
        roots: list[float] = []
        for i in range(len(g) - 1):
            if not (np.isfinite(d[i]) and np.isfinite(d[i + 1])):
                continue
            if d[i] == 0.0:
                r = float(g[i])
            elif d[i] * d[i + 1] < 0:
                try:
                    r = float(brentq(lambda x: float(fa(x) - fb(x)), g[i], g[i + 1]))
                except Exception:
                    continue
            else:
                continue
            if not roots or abs(r - roots[-1]) > 1e-8:
                roots.append(r)
        if not roots:
            out.append((int(a), int(b), float("nan"), 0, ""))
            continue
        chosen = min(roots, key=lambda r: abs(r - float(tc)))
        out.append((int(a), int(b), float(chosen), int(len(roots)),
                    ";".join(f"{r:.6f}" for r in roots)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--spec", type=Path, required=True)
    ap.add_argument("--package-src", type=Path, required=True)
    ap.add_argument("--tc-table", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ns = ap.parse_args()

    print(f"FIGURE_BUILD_VERSION = {FIGURE_BUILD_VERSION}")
    print(f"script sha256        = {_sha(Path(__file__).resolve())}")

    spec_path = ns.spec.resolve()
    if not spec_path.exists():
        raise SystemExit(f"FAIL-CLOSED: --spec does not exist: {spec_path}")
    try:
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise SystemExit(f"FAIL-CLOSED: --spec is not readable JSON: {exc}")
    if spec.get("spec_name") != "FGT_CORRECTION_AWARE_CRITICAL_AUDIT_v3.2.1":
        raise SystemExit("FAIL-CLOSED: --spec is not the locked v3.2.1 SPEC_LOCK.json")

    src = ns.package_src.resolve()
    if not (src / "fgt_csaudit").is_dir():
        raise SystemExit(f"FAIL-CLOSED: --package-src has no fgt_csaudit: {src}")
    sys.path.insert(0, str(src))
    from fgt_csaudit.io import validate_project_inputs, load_all_cases  # noqa: E402

    if ns.tc_table is not None:
        tc_path = ns.tc_table.resolve()
        if not tc_path.exists():
            raise SystemExit(f"FAIL-CLOSED: --tc-table does not exist: {tc_path}")
    else:
        hits = sorted(ns.results.resolve().rglob("TABLE_TC_RG_CROSSING_BOOTSTRAP.csv"))
        if not hits:
            raise SystemExit("FAIL-CLOSED: TABLE_TC_RG_CROSSING_BOOTSTRAP.csv not found")
        if len(hits) > 1:
            listing = "\n  ".join(str(h) for h in hits)
            raise SystemExit("FAIL-CLOSED: more than one TABLE_TC_RG_CROSSING_BOOTSTRAP.csv found; "
                             f"pass --tc-table explicitly.\n  {listing}")
        tc_path = hits[0]
    tc_tab = pd.read_csv(tc_path)
    locked_tc = {str(r.case_label): float(r.tc_bootstrap_median) for r in tc_tab.itertuples()}

    paths, _ = validate_project_inputs(ns.project_root.resolve(), spec)
    cases = load_all_cases(paths, spec)

    diag: list[dict] = []
    fig, axes = plt.subplots(2, 2, figsize=(6.9, 5.0))
    for i, (label, title, half) in enumerate(ROWS):
        if label not in locked_tc:
            raise SystemExit(f"FAIL-CLOSED: no locked Tc row for {label}")
        case = cases[label]
        cur = case.central_curves()
        tc = locked_tc[label]
        cmap = plt.cm.viridis(np.linspace(0.12, 0.88, len(case.sizes)))
        for j, (ch, ylab) in enumerate(COLS):
            ax = axes[i, j]
            stack = []
            for col, L in zip(cmap, case.sizes):
                t, y = cur[ch][L]
                m = (t >= tc - half) & (t <= tc + half)
                ax.plot(t[m], y[m], "-o", ms=3.2, lw=1.1, color=col,
                        mfc="white", mew=0.8, label=f"$L={L}$")
                stack.append(y[m])
            ax.axvline(tc, color="#c00000", lw=1.0, ls="--", zorder=1)
            ax.set_xlim(tc - half, tc + half)
            v = np.concatenate(stack)
            lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
            ax.set_ylim(lo - 0.04 * (hi - lo), hi + 0.26 * (hi - lo))
            ax.set_ylabel(ylab)
            ax.grid(alpha=0.22)
            if i == len(ROWS) - 1:
                ax.set_xlabel(r"Temperature, $T/(J/k_B)$")
            if j == 0:
                ax.text(0.03, 0.96, title, transform=ax.transAxes, va="top", fontsize=8.5)
            for a, b, root, nroot, allroots in _adjacent_crossings(
                    cur[ch], case.sizes, tc - half, tc + half, tc):
                diag.append({"case_label": label, "channel": ch, "L_small": a, "L_large": b,
                             "crossing_temperature": root, "n_roots_in_window": nroot,
                             "all_roots_in_window": allroots, "locked_tc": tc,
                             "window_half_width": half,
                             "root_selection_rule": "closest_to_locked_tc"})
    axes[0, 1].legend(loc="best", ncol=2, framealpha=0.94, edgecolor="0.55",
                      fancybox=False, handlelength=1.4, labelspacing=0.28, borderpad=0.38)
    fig.tight_layout()

    out = ns.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cat = []
    for fmt in FORMATS:
        p = out / f"{STEM}.{fmt}"
        fig.savefig(p, **({"dpi": DPI} if fmt == "png" else {}))
        cat.append({"figure": STEM, "format": fmt, "path": str(p),
                    "sha256": _sha(p), "bytes": p.stat().st_size})
    plt.close(fig)

    d = pd.DataFrame(diag)
    summary = []
    for (lab, ch), g in d.groupby(["case_label", "channel"]):
        v = g.crossing_temperature.to_numpy(float)
        v = v[np.isfinite(v)]
        nmax = int(g.n_roots_in_window.max()) if len(g) else 0
        summary.append({
            "case_label": lab, "channel": ch,
            "n_adjacent_pairs_attempted": int(len(g)),
            "n_adjacent_pairs_with_crossing": int(len(v)),
            "max_roots_in_any_pair": nmax,
            "any_pair_multi_root": bool(nmax > 1),
            "crossing_min": float(v.min()) if len(v) else float("nan"),
            "crossing_max": float(v.max()) if len(v) else float("nan"),
            "crossing_spread": float(v.max() - v.min()) if len(v) else float("nan")})
    sdf = pd.DataFrame(summary)

    def _pairset(sub: pd.DataFrame) -> str:
        ok = sub.loc[np.isfinite(sub.crossing_temperature.to_numpy(float))]
        return ";".join(f"{int(a)}-{int(b)}" for a, b in
                        sorted(zip(ok.L_small, ok.L_large)))

    pairsets = {(str(lab), str(ch)): _pairset(g)
                for (lab, ch), g in d.groupby(["case_label", "channel"])}
    sdf["contributing_pair_set"] = [pairsets[(r.case_label, r.channel)]
                                    for r in sdf.itertuples()]
    by_ch = sdf.groupby("channel")["contributing_pair_set"].nunique()
    sdf["comparable_across_cases"] = sdf["channel"].map(lambda c: bool(by_ch[c] == 1))
    sdf["comparability_note"] = sdf["comparable_across_cases"].map(
        {True: "identical set of contributing adjacent size pairs in all cases",
         False: "contributing adjacent size pairs differ across cases; do not quote a ratio"})
    d.to_csv(out / "TABLE_CROSSING_DRIFT_DIAGNOSTIC.csv", index=False)
    sdf.to_csv(out / "TABLE_CROSSING_DRIFT_SUMMARY.csv", index=False)
    summary = sdf.to_dict(orient="records")

    (out / f"{STEM}.metadata.json").write_text(json.dumps({
        "figure": STEM,
        "figure_build_version": FIGURE_BUILD_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "spec_source": str(spec_path), "spec_sha256": _sha(spec_path),
        "locked_tc_source": str(tc_path), "locked_tc_sha256": _sha(tc_path),
        "locked_tc_used": locked_tc,
        "legacy_estimator_used": False,
        "tc_recomputed_in_this_script": False,
        "root_selection_rule": "closest_to_locked_tc (fixed before inspection)",
        "crossing_counts_reported_per_pair": True,
        "comparability_rule": "identical contributing adjacent-pair set across cases",
        "plotted_window_half_widths": {lab: half for lab, _, half in ROWS},
        "plotted_quantity": "central disorder-averaged curves (CaseData.central_curves)",
        "uncertainty_shown_in_figure": False,
        "files": cat,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    with (out / f"{STEM}_CATALOG.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["figure", "format", "path", "sha256", "bytes"])
        w.writeheader(); w.writerows(cat)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
