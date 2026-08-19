from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .provenance import atomic_write_json, sha256_file


def _save(fig: Any, fig_dir: Path, stem: str, spec: dict[str, Any], metadata: dict[str, Any]) -> list[dict[str, Any]]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for fmt in spec["figures"]["formats"]:
        p = fig_dir / f"{stem}.{fmt}"
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = int(spec["figures"]["png_dpi"])
        fig.savefig(p, **kwargs)
        rows.append({"figure": stem, "format": fmt, "path": str(p), "sha256": sha256_file(p), "bytes": p.stat().st_size})
    atomic_write_json(fig_dir / f"{stem}.metadata.json", metadata)
    plt.close(fig)
    return rows


def build_figures(summary: dict[str, pd.DataFrame], decisions: pd.DataFrame, fig_dir: Path, spec: dict[str, Any]) -> pd.DataFrame:
    plt.rcParams.update({"font.size": float(spec["figures"]["font_size"]), "legend.fontsize": float(spec["figures"]["legend_font_size"])})
    cat: list[dict[str, Any]] = []

    tc = summary["tc"].sort_values("p")
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    p = tc.p.to_numpy(float); y = tc.tc_bootstrap_median.to_numpy(float)
    elo = y - tc.tc_bootstrap_ci_low.to_numpy(float); ehi = tc.tc_bootstrap_ci_high.to_numpy(float) - y
    ax.errorbar(p, y, yerr=np.vstack([elo, ehi]), fmt="o-", capsize=3, label="Joint RG-invariant crossing Tc; 95% quenched bootstrap CI")
    leg = [spec["legacy_tc_provenance"][x]["tc"] for x in tc.case_label]
    ax.plot(p, leg, "x--", label="Legacy susceptibility-shift Tc (nu=1 conditioned; provenance only)")
    p1 = tc.loc[np.isclose(tc.p, 1.0)]
    if len(p1):
        ax.scatter([1.0], [float(spec["exact_honeycomb_ising_tc"])], marker="*", s=90, label="Exact pristine honeycomb Ising Tc")
    ax.set_xlabel("Site occupation probability, p")
    ax.set_ylabel(r"Critical temperature, $T_c/(J/k_B)$")
    ax.set_title("Tc from dimensionless RG-invariant crossings")
    ax.legend()
    ax.grid(alpha=0.25)
    cat += _save(fig, fig_dir, "FIG01_rg_crossing_tc_vs_p", spec, {"meaning": "Joint Binder and xi/L crossing Tc with bootstrap CI; legacy fixed-nu Tc shown only for provenance."})

    nu = summary["nu"].sort_values(["p", "channel"])
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    offsets = {"binder_roa": -0.006, "xi_over_L": 0.006}
    labels = {"binder_roa": "Binder ratio-of-averages", "xi_over_L": r"Correlation-length ratio $\xi/L$"}
    for ch in ("binder_roa", "xi_over_L"):
        g = nu.loc[nu.channel == ch].sort_values("p")
        x = g.p.to_numpy(float) + offsets[ch]
        med = g.nu_bootstrap_median.to_numpy(float)
        lo = g.nu_bootstrap_ci_low.to_numpy(float); hi = g.nu_bootstrap_ci_high.to_numpy(float)
        ax.errorbar(x, med, yerr=np.vstack([med-lo, hi-med]), fmt="o", capsize=3, label=labels[ch] + "; bootstrap CI")
        for xx, rrlo, rrhi in zip(x, g.nu_robustness_low, g.nu_robustness_high):
            ax.plot([xx, xx], [rrlo, rrhi], linewidth=3, alpha=0.35)
    ax.axhline(1.0, linestyle="--", linewidth=1, label=r"2D Ising reference $\nu=1$")
    ax.set_xlabel("Site occupation probability, p")
    ax.set_ylabel(r"Correlation-length exponent estimate, $\nu$")
    ax.set_title("Statistical intervals and correction-aware robustness envelopes")
    ax.legend()
    ax.grid(alpha=0.25)
    cat += _save(fig, fig_dir, "FIG02_nu_bootstrap_and_robustness", spec, {"meaning": "Thin error bars are quenched-realization percentile bootstrap intervals. Thick translucent spans are conservative robustness envelopes and are not confidence intervals."})

    dr = summary["drift"].sort_values(["p", "channel"])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for ch in ("binder_roa", "xi_over_L"):
        g = dr.loc[dr.channel == ch].sort_values("p")
        x = g.p.to_numpy(float) + offsets[ch]
        med = g.delta_nu_drop_minus_full_median.to_numpy(float)
        lo = g.delta_ci_low.to_numpy(float); hi = g.delta_ci_high.to_numpy(float)
        ax.errorbar(x, med, yerr=np.vstack([med-lo, hi-med]), fmt="o", capsize=3, label=labels[ch])
    ax.axhline(0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("Site occupation probability, p")
    ax.set_ylabel(r"Paired drift, $\Delta\nu=\nu_{L\geq60}-\nu_{L\geq40}$")
    ax.set_title("Paired finite-size drift diagnostic")
    ax.legend()
    ax.grid(alpha=0.25)
    cat += _save(fig, fig_dir, "FIG03_paired_Lmin_drift", spec, {"meaning": "Paired bootstrap drift from the same disorder resample; drift toward nu=1 is a veto against overinterpreting effective exponents."})

    ra = summary["ratios"].sort_values(["p", "channel"])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ratio_labels = {"abs_m": r"$\beta/\nu$ consistency ratio", "chi_abs": r"$\gamma/\nu$ consistency ratio"}
    for ch in ("abs_m", "chi_abs"):
        g = ra.loc[ra.channel == ch].sort_values("p")
        x = g.p.to_numpy(float) + offsets.get("binder_roa" if ch == "abs_m" else "xi_over_L", 0.0)
        ref = g.ising_reference.to_numpy(float)
        med = g.ratio_median.to_numpy(float) / ref
        lo = g.ratio_ci_low.to_numpy(float) / ref; hi = g.ratio_ci_high.to_numpy(float) / ref
        ax.errorbar(x, med, yerr=np.vstack([med-lo, hi-med]), fmt="o", capsize=3, label=ratio_labels[ch])
    ax.axhline(1.0, linestyle="--", linewidth=1, label="2D Ising reference")
    ax.set_xlabel("Site occupation probability, p")
    ax.set_ylabel("Estimated exponent ratio / 2D Ising reference")
    ax.set_title("Exponent-ratio consistency at independently estimated Tc")
    ax.legend()
    ax.grid(alpha=0.25)
    cat += _save(fig, fig_dir, "FIG04_exponent_ratio_consistency", spec, {"meaning": "Diagnostic only; exponent-ratio agreement is not used to tune Tc and is not proof of universality."})

    dec = decisions.sort_values("p")
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ymap = {"NU1_COMPATIBLE": 3, "INCONCLUSIVE_CORRECTION_DOMINATED": 2, "INCONCLUSIVE_LIMITED_RANGE": 1, "EVIDENCE_AGAINST_NU1": 0}
    vals = [ymap.get(x, -1) for x in dec.decision]
    ax.scatter(dec.p, vals, s=55)
    ax.set_yticks([0,1,2,3], ["Evidence against nu=1", "Inconclusive: limited range", "Inconclusive: corrections", "nu=1 compatible"])
    ax.set_xlabel("Site occupation probability, p")
    ax.set_title("Fail-closed critical-scaling decision")
    ax.grid(axis="x", alpha=0.25)
    cat += _save(fig, fig_dir, "FIG05_decision_summary", spec, {"meaning": "Decision categories are framework-specific and never constitute proof of a new universality class."})

    return pd.DataFrame(cat)
