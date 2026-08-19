#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BUILD_PUBLICATION_FIGURES_v321_1.py
===================================
يعيد بناء الأشكال من الجداول المصحَّحة. قراءة فقط للمدخلات؛ لا bootstrap،
لا إعادة تشغيل 03، ولا مساس بمجلد نتائج v3.2.1 الأصلي.

الاستعمال:
    python BUILD_PUBLICATION_FIGURES_v321_1.py ^
        --original "<...>\\final_csaudit_v3_2_1_correction_aware" ^
        --corrected "<...>\\final_csaudit_v3_2_1_1_decision_corrected"

يُنتج داخل مجلد المصحَّح: figures_v321_1\\  بصيغ PDF (متجهي) و PNG 600dpi و SVG،
مع FIGURE_CATALOG_v321_1.csv يحمل SHA-256 لكل ملف.

ما الذي تغيّر عن أشكال v3.2.1، ولماذا
-------------------------------------
FIG01  أُضيفت لوحة بواقٍ سفلية.  في الشكل الأصلي يمتد المحور من 0.9 إلى 1.5
       فتصبح إزاحة Tc البالغة 0.003 والفروق عن القيمة المضبوطة غير مرئية
       بصرياً.  لوحة البواقي تُظهر الكمية التي تحمل المعنى فعلاً.
FIG02  أُصلح خطآن: (أ) أغلفة المتانة كانت تُرسم بألوان دورية عشوائية لا تطابق
       ألوان القنوات ولا يوجد لها مدخل في المفتاح؛ (ب) لم تكن تُظهر نسبة
       الارتطام، وهي الكمية الحاكمة للقرار بعد التصحيح.  والغلاف يُرسم الآن
       بحدّ متقطّع لا كتلة مصمتة، لأنه ليس فترة ثقة.
FIG05  أُعيد توليده بالكامل من TABLE_FINAL_DECISIONS_v321_1.csv.  الشكل الأصلي
       يعرض NU1_COMPATIBLE للحالات الأربع، وثلاثة منها لم تعد قائمة.
FIG03  محتفظ به مع تحسينات عرض وتعليم واضح للحالات المُستبعَدة ببوابة.
FIG04  محتفظ به مع تعليم القيمة المرجعية.
FIG06  شكل جديد: أربعة تشخيصات متكاملة (complementary) لقابلية تحديد nu.  هذا هو الشكل الذي
       يحمل النتيجة المنهجية للورقة.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

FIGURE_BUILD_VERSION = "v321_1-figs-2026-08-12-legend-inside-uniform"

EXACT_TC = 1.518651435000414
CH = ("binder_roa", "xi_over_L")
CHLAB = {"binder_roa": "Binder ratio-of-averages",
         "xi_over_L": r"Correlation-length ratio $\xi/L$"}
CHCOL = {"binder_roa": "#1f4e79", "xi_over_L": "#a5390d"}
OFF = {"binder_roa": -0.007, "xi_over_L": 0.007}
FORMATS = ("pdf", "png", "svg")
DPI = 600

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.8, "grid.linewidth": 0.5, "lines.linewidth": 1.2,
    "figure.dpi": 120, "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})


LEGEND_KW = dict(framealpha=0.94, fancybox=False, edgecolor="0.55",
                 handlelength=1.6, handletextpad=0.55, labelspacing=0.36,
                 borderpad=0.42, borderaxespad=0.5)
LEGEND_FS = 7.0
LEGEND_FS_SMALL = 6.5


def _place_legend(ax, *, fontsize=LEGEND_FS, headroom=0.0, handles=None,
                  labels=None, prefer=("upper left", "upper right",
                                       "lower left", "lower right")):
    """Put a compact legend inside the axes in the emptiest corner.

    ``headroom`` reserves a fraction of the current y-span above the data so the
    legend box has somewhere to sit that cannot collide with a marker or cap.
    matplotlib's 'best' placement then picks the free corner; ``prefer`` is used
    only when 'best' would be ambiguous.
    """
    if headroom:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + headroom * (hi - lo))
    kw = dict(LEGEND_KW, fontsize=fontsize)
    if handles is not None:
        return ax.legend(handles=handles, labels=labels, loc="best", **kw)
    return ax.legend(loc="best", **kw)


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _find(root: Path, name: str) -> Path:
    for c in (root / name, root / "tables" / name, root / "manifests" / name):
        if c.exists():
            return c
    hits = sorted(root.rglob(name))
    if not hits:
        raise FileNotFoundError(f"{name} not found under {root}")
    return hits[0]


def _save(fig, out: Path, stem: str, meaning: str, cat: list[dict]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for fmt in FORMATS:
        p = out / f"{stem}.{fmt}"
        fig.savefig(p, **({"dpi": DPI} if fmt == "png" else {}))
        cat.append({"figure": stem, "format": fmt, "path": str(p),
                    "sha256": _sha(p), "bytes": p.stat().st_size})
    (out / f"{stem}.metadata.json").write_text(
        json.dumps({"figure": stem, "meaning": meaning,
                    "figure_build_version": FIGURE_BUILD_VERSION,
                    "decision_layer_version": "v3.2.1.1-decision-corrected"},
                   indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig01(tc: pd.DataFrame, legacy: dict[str, float], out: Path, cat: list) -> None:
    """Tc(p) with a residual panel that makes the sub-percent structure visible."""
    d = tc.sort_values("p")
    p = d.p.to_numpy(float)
    y = d.tc_bootstrap_median.to_numpy(float)
    lo = d.tc_bootstrap_ci_low.to_numpy(float)
    hi = d.tc_bootstrap_ci_high.to_numpy(float)
    leg = np.array([legacy[c] for c in d.case_label], float)

    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 4.3), sharex=True,
                               gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.08})

    a.errorbar(p, y, yerr=np.vstack([y - lo, hi - y]), fmt="o-", ms=4.5, capsize=2.5,
               color=CHCOL["binder_roa"], mfc="white", mew=1.2, zorder=3,
               label="Joint RG-invariant crossing $T_c$\n(95% quenched bootstrap CI)")
    a.plot(1.0, EXACT_TC, marker="*", ms=11, color="#c00000", ls="none", zorder=4,
           label="Exact honeycomb Ising $T_c$ ($p=1$)")
    a.set_ylabel(r"$T_c\;/\;(J/k_B)$")
    a.grid(alpha=0.22)
    _place_legend(a, headroom=0.30)

    b.axhline(0.0, color="0.35", lw=0.9, ls="--", zorder=1)
    b.errorbar(p, (y - leg) * 1e3, yerr=np.vstack([y - lo, hi - y]) * 1e3, fmt="s",
               ms=4, capsize=2.5, color="#2a6f4e", mfc="white", mew=1.1, zorder=3,
               label=r"$T_c^{\rm RG}-T_c^{\rm legacy}$")
    b.plot(1.0, (EXACT_TC - legacy["pristine_p100"]) * 1e3, marker="*", ms=11,
           color="#c00000", ls="none", zorder=4,
           label=r"$T_c^{\rm exact}-T_c^{\rm legacy}$")
    b.set_xlabel("Site occupation probability, $p$")
    b.set_ylabel(r"$\Delta T_c\;\times10^{3}$")
    b.grid(alpha=0.22)
    b.margins(y=0.55)
    _place_legend(b)

    _save(fig, out, "FIG01_rg_crossing_tc_vs_p", (
        "Panel (a): joint Binder and xi/L crossing Tc with 95% quenched-realization "
        "percentile bootstrap intervals; the exact honeycomb-Ising value is shown at p=1 "
        "as an external validation point. Panel (b): residual against the legacy "
        "fixed-nu susceptibility-shift Tc, which is provenance only. The residual panel "
        "is required because the sub-percent structure is not resolvable on the absolute "
        "scale of panel (a)."), cat)


def fig02(nu: pd.DataFrame, out: Path, cat: list) -> None:
    """nu intervals, envelopes drawn as non-CI outlines, and the railing fraction."""
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.6, 4.7), sharex=True,
                               gridspec_kw={"height_ratios": [2.5, 1.0], "hspace": 0.10})
    a.axhline(1.0, ls="--", lw=1.0, color="0.30", zorder=1,
              label=r"2D Ising reference $\nu=1$")
    for ch in CH:
        g = nu.loc[nu.channel == ch].sort_values("p")
        x = g.p.to_numpy(float) + OFF[ch]
        med = g.nu_bootstrap_median.to_numpy(float)
        lo = g.nu_bootstrap_ci_low.to_numpy(float)
        hi = g.nu_bootstrap_ci_high.to_numpy(float)
        rl = g.nu_robustness_low.to_numpy(float)
        rh = g.nu_robustness_high.to_numpy(float)
        # Envelope: dashed outline in the channel colour. Not a confidence interval,
        # so it is deliberately not drawn as a solid filled span.
        for xx, l, h in zip(x, rl, rh):
            a.plot([xx, xx], [l, h], color=CHCOL[ch], lw=0.8, ls=(0, (2, 1.6)),
                   alpha=0.75, zorder=2)
            for yy in (l, h):
                a.plot([xx - 0.006, xx + 0.006], [yy, yy], color=CHCOL[ch],
                       lw=0.8, alpha=0.75, zorder=2)
        a.errorbar(x, med, yerr=np.vstack([med - lo, hi - med]), fmt="o", ms=4.5,
                   capsize=2.5, color=CHCOL[ch], mfc="white", mew=1.2, zorder=4,
                   label=CHLAB[ch] + " (bootstrap CI)")
    a.set_ylabel(r"Correlation-length exponent, $\nu$")
    a.grid(alpha=0.22)
    h, l = a.get_legend_handles_labels()
    h.append(Line2D([], [], color="0.35", lw=0.9, ls=(0, (2, 1.6))))
    l.append("Robustness envelope (NOT a CI)")
    # Reserve headroom so the four-entry legend remains inside the axes without
    # colliding with the tallest robustness envelope.
    ymax = max(float(nu.nu_robustness_high.max()), float(nu.nu_bootstrap_ci_high.max()))
    ymin = min(float(nu.nu_robustness_low.min()), float(nu.nu_bootstrap_ci_low.min()))
    span = ymax - ymin
    a.set_ylim(ymin - 0.04 * span, ymax + 0.40 * span)
    _place_legend(a, fontsize=LEGEND_FS_SMALL, handles=h, labels=l)

    thr = 0.10
    for ch in CH:
        g = nu.loc[nu.channel == ch].sort_values("p")
        b.plot(g.p.to_numpy(float) + OFF[ch],
               g.bootstrap_boundary_fraction.to_numpy(float), "o-", ms=4.5,
               color=CHCOL[ch], mfc="white", mew=1.2, label=CHLAB[ch])
    b.axhline(thr, color="#c00000", lw=1.0, ls="--")
    b.text(0.802, thr + 0.02, "locked gate = 0.10", color="#c00000", fontsize=7)
    b.set_xlabel("Site occupation probability, $p$")
    b.set_ylabel("Boundary-hit\nfraction")
    b.set_ylim(-0.03, 0.56)
    b.grid(alpha=0.22)
    _place_legend(b)

    _save(fig, out, "FIG02_nu_intervals_and_identifiability", (
        "Panel (a): quenched-realization percentile bootstrap intervals for nu (solid "
        "markers and caps) and the conservative correction-aware robustness envelope "
        "(dashed outline). The envelope is not a confidence interval and is used only to "
        "make EVIDENCE_AGAINST_NU1 harder to reach. Panel (b): fraction of bootstrap "
        "draws whose nu estimate rails at the symmetric feasible search bound. Values "
        "above the locked 0.10 gate mean nu is not identified and the case is labelled "
        "INCONCLUSIVE_LIMITED_RANGE."), cat)


def fig03(dr: pd.DataFrame, out: Path, cat: list) -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.7))
    ax.axhline(0.0, ls="--", lw=1.0, color="0.30", zorder=1)
    for ch in CH:
        g = dr.loc[dr.channel == ch].sort_values("p")
        x = g.p.to_numpy(float) + OFF[ch]
        m = g.delta_nu_drop_minus_full_median.to_numpy(float)
        lo = g.delta_ci_low.to_numpy(float)
        hi = g.delta_ci_high.to_numpy(float)
        ax.errorbar(x, m, yerr=np.vstack([m - lo, hi - m]), fmt="o", ms=4.5, capsize=2.5,
                    color=CHCOL[ch], mfc="white", mew=1.2, label=CHLAB[ch])
    ax.set_xlabel("Site occupation probability, $p$")
    ax.set_ylabel(r"$\Delta\nu=\nu_{L\geq60}-\nu_{L\geq40}$")
    ax.grid(alpha=0.22)
    _place_legend(ax, headroom=0.26)
    _save(fig, out, "FIG03_paired_Lmin_drift", (
        "Paired finite-size drift computed from the same disorder resample in each "
        "bootstrap draw. No interval excludes zero, so the correction-dominated veto is "
        "not triggered for any case; the diluted cases are gated out by the "
        "identifiability criterion instead."), cat)


def fig04(ra: pd.DataFrame, out: Path, cat: list) -> None:
    lab = {"abs_m": r"$\beta/\nu$", "chi_abs": r"$\gamma/\nu$"}
    col = {"abs_m": CHCOL["binder_roa"], "chi_abs": CHCOL["xi_over_L"]}
    off = {"abs_m": -0.007, "chi_abs": 0.007}
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    ax.axhline(1.0, ls="--", lw=1.0, color="0.30", zorder=1, label="2D Ising reference")
    for ch in ("abs_m", "chi_abs"):
        g = ra.loc[ra.channel == ch].sort_values("p")
        ref = g.ising_reference.to_numpy(float)
        x = g.p.to_numpy(float) + off[ch]
        m = g.ratio_median.to_numpy(float) / ref
        lo = g.ratio_ci_low.to_numpy(float) / ref
        hi = g.ratio_ci_high.to_numpy(float) / ref
        ax.errorbar(x, m, yerr=np.vstack([m - lo, hi - m]), fmt="o", ms=4.5, capsize=2.5,
                    color=col[ch], mfc="white", mew=1.2, label=lab[ch] + " / reference")
    ax.set_xlabel("Site occupation probability, $p$")
    ax.set_ylabel("Estimated ratio / 2D Ising value")
    ax.grid(alpha=0.22)
    lo_all, hi_all = [], []
    for ch in ("abs_m", "chi_abs"):
        g = ra.loc[ra.channel == ch]
        ref = g.ising_reference.to_numpy(float)
        lo_all += list(g.ratio_ci_low.to_numpy(float) / ref)
        hi_all += list(g.ratio_ci_high.to_numpy(float) / ref)
    ymin, ymax = min(lo_all), max(hi_all)
    span = ymax - ymin
    ax.set_ylim(ymin - 0.05 * span, ymax + 0.30 * span)
    _place_legend(ax)
    _save(fig, out, "FIG04_exponent_ratio_consistency", (
        "Exponent ratios evaluated at the independently estimated RG-crossing Tc, "
        "normalised by their exact 2D Ising values. Consistency diagnostic only: these "
        "ratios are never used to tune Tc and do not by themselves establish "
        "universality."), cat)


def fig05(dec: pd.DataFrame, out: Path, cat: list) -> None:
    order = ["EVIDENCE_AGAINST_NU1", "INCONCLUSIVE_LIMITED_RANGE",
             "INCONCLUSIVE_CORRECTION_DOMINATED", "NU1_COMPATIBLE"]
    pretty = {"EVIDENCE_AGAINST_NU1": "Evidence against $\\nu=1$",
              "INCONCLUSIVE_LIMITED_RANGE": "Inconclusive:\n$\\nu$ not identified",
              "INCONCLUSIVE_CORRECTION_DOMINATED": "Inconclusive:\ncorrection dominated",
              "NU1_COMPATIBLE": "$\\nu=1$ compatible"}
    ymap = {k: i for i, k in enumerate(order)}
    d = dec.sort_values("p")
    fig, ax = plt.subplots(figsize=(3.9, 2.6))
    for _, r in d.iterrows():
        yy = ymap.get(r.decision, -1)
        gated = not bool(r.primary_boundary_gate_pass)
        ax.scatter([r.p], [yy], s=70, zorder=3,
                   color="#c00000" if gated else "#2a6f4e",
                   marker="X" if gated else "o", edgecolor="white", linewidth=0.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([pretty[k] for k in order])
    ax.set_ylim(-0.55, len(order) - 0.05)
    ax.set_xlim(0.78, 1.02)
    ax.set_xlabel("Site occupation probability, $p$")
    ax.grid(axis="x", alpha=0.22)
    _place_legend(ax, handles=[
        Line2D([], [], ls="none", marker="o", ms=6, color="#2a6f4e"),
        Line2D([], [], ls="none", marker="X", ms=7, color="#c00000")],
        labels=["identifiability gate passed", "gated out (boundary railing)"])
    _save(fig, out, "FIG05_decision_summary", (
        "Fail-closed decision per case under the corrected v3.2.1.1 decision layer. "
        "Categories are framework-specific and never constitute proof of a new "
        "universality class. Crossed markers denote cases removed by the locked "
        "boundary-railing gate, for which no nu statement is supported in either "
        "direction."), cat)


def fig06(dec: pd.DataFrame, nu: pd.DataFrame, out: Path, cat: list) -> None:
    """Four complementary identifiability diagnostics: the methodological result."""
    d = dec.sort_values("p")
    p = d.p.to_numpy(float)
    rail = np.array([max(float(r.binder_bootstrap_boundary_fraction),
                         float(r.xi_bootstrap_boundary_fraction)) for _, r in d.iterrows()])
    dnu = d.channel_delta_nu_ci_width.to_numpy(float)
    frac = np.array([max(float(r.binder_nu_ci_width_over_search_width),
                         float(r.xi_nu_ci_width_over_search_width)) for _, r in d.iterrows()])
    width = np.array([max(float(r.binder_nu_ci_width), float(r.xi_nu_ci_width))
                      for _, r in d.iterrows()])

    panels = [
        (rail, "Boundary-hit fraction\n(worst channel)", 0.10, "locked gate"),
        (dnu, r"Channel disagreement" "\n" r"width of $\Delta\nu$ CI", None, None),
        (frac, "Bootstrap CI width /\nfeasible search width", None, None),
        (width, r"Bootstrap CI width in $\nu$" "\n(worst channel)", None, None),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 4.2))
    for ax, (v, lab, thr, thrlab) in zip(axes.ravel(), panels):
        ax.plot(p, v, "o-", ms=5, color=CHCOL["binder_roa"], mfc="white", mew=1.2, zorder=3)
        ax.scatter([1.0], [v[-1]], s=90, marker="*", color="#2a6f4e", zorder=4)
        if thr is not None:
            ax.axhline(thr, color="#c00000", lw=1.0, ls="--")
            ax.text(0.80, thr * 1.12, thrlab, color="#c00000", fontsize=7)
        ax.set_ylabel(lab, fontsize=8)
        ax.set_xlabel("$p$", fontsize=8)
        ax.grid(alpha=0.22)
        ax.set_xlim(0.78, 1.03)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.16 * (hi - lo))
    axes.ravel()[0].legend(handles=[
        Line2D([], [], ls="-", marker="o", ms=5, color=CHCOL["binder_roa"],
               mfc="white", mew=1.2),
        Line2D([], [], ls="none", marker="*", ms=9, color="#2a6f4e")],
        labels=["diluted cases", "pristine $p=1$ (calibrated)"],
        loc="best", **dict(LEGEND_KW, fontsize=6.6))
    fig.suptitle("Four complementary identifiability diagnostics for $\\nu$", fontsize=10, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    _save(fig, out, "FIG06_identifiability_metrics", (
        "Four diagnostics of distinct construction, all computed from the same production "
        "bootstrap and therefore complementary rather than statistically independent. All "
        "four separate the analytically calibrated pristine case (green star) from the "
        "diluted cases by a large margin, indicating that the limitation is the accessible "
        "size range rather than the estimation method."), cat)


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--corrected", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--spec", type=Path, required=True,
                    help="Path to the SPEC_LOCK.json actually used by the audit run. Mandatory.")
    ns = ap.parse_args()

    print(f"FIGURE_BUILD_VERSION = {FIGURE_BUILD_VERSION}")
    print(f"script file          = {Path(__file__).resolve()}")
    print(f"script sha256        = {_sha(Path(__file__).resolve())}")
    orig, corr = ns.original.resolve(), ns.corrected.resolve()
    out = (ns.out or corr / "figures_v321_1").resolve()

    tc = pd.read_csv(_find(orig, "TABLE_TC_RG_CROSSING_BOOTSTRAP.csv"))
    nu = pd.read_csv(_find(orig, "TABLE_NU_CORRECTION_AWARE.csv"))
    dr = pd.read_csv(_find(orig, "TABLE_PAIRED_LMIN_DRIFT.csv"))
    ra = pd.read_csv(_find(orig, "TABLE_EXPONENT_RATIO_CONSISTENCY.csv"))
    dec = pd.read_csv(_find(corr, "TABLE_FINAL_DECISIONS_v321_1.csv"))

    spec_path = ns.spec.resolve()
    if not spec_path.exists():
        raise SystemExit(f"FAIL-CLOSED: --spec path does not exist: {spec_path}")
    try:
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise SystemExit(f"FAIL-CLOSED: --spec is not readable JSON: {spec_path}: {exc}")
    if "legacy_tc_provenance" not in spec:
        raise SystemExit(f"FAIL-CLOSED: --spec lacks legacy_tc_provenance: {spec_path}")
    legacy = {k: float(v["tc"]) for k, v in spec["legacy_tc_provenance"].items()}
    spec_record = {"spec_source": str(spec_path), "spec_sha256": _sha(spec_path)}
    print(f"spec provenance: {spec_record['spec_source']}")

    cat: list[dict] = []
    fig01(tc, legacy, out, cat)
    fig02(nu, out, cat)
    fig03(dr, out, cat)
    fig04(ra, out, cat)
    fig05(dec, out, cat)
    fig06(dec, nu, out, cat)

    (out / "FIGURE_PROVENANCE.json").write_text(
        json.dumps({"tool": "BUILD_PUBLICATION_FIGURES_v321_1",
                    "figure_build_version": FIGURE_BUILD_VERSION,
                    "script_path": str(Path(__file__).resolve()),
                    "script_sha256": _sha(Path(__file__).resolve()),
                    "decision_layer_version": "v3.2.1.1-decision-corrected",
                    "original_results_dir": str(orig), "corrected_results_dir": str(corr),
                    **spec_record,
                    "figures": sorted({c["figure"] for c in cat})},
                   indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    cpath = out / "FIGURE_CATALOG_v321_1.csv"
    with cpath.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["figure", "format", "path", "sha256", "bytes"])
        w.writeheader()
        w.writerows(cat)

    print(f"\nwrote {len(cat)} figure files ({len(cat)//len(FORMATS)} figures) to:\n  {out}")
    for stem in sorted({c['figure'] for c in cat}):
        print(f"  - {stem}")
    print(f"\ncatalog with SHA-256: {cpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
