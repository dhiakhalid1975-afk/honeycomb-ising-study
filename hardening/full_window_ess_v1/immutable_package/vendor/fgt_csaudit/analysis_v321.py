from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .correction_fit import (
    exponent_ratio_at_tc,
    fit_nu_fixed_tc,
    restrict_curves,
    symmetric_feasible_nu_bounds,
)
from .pb import build_locked_support
from .rg_tc import estimate_joint_tc

PRIMARY = ("binder_roa", "xi_over_L")


def _window_curves(curves_all: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]], sizes: list[int] | tuple[int, ...]) -> dict[str, dict[int, tuple[np.ndarray, np.ndarray]]]:
    return {ch: restrict_curves(curves_all[ch], sizes) for ch in curves_all}


def build_central_case(case: Any, spec: dict[str, Any]) -> dict[str, Any]:
    curves_all = case.central_curves()
    legacy = spec["legacy_tc_provenance"][case.label]
    branch_center = float(legacy["tc"])
    windows = {
        "full": list(spec["size_windows"]["full"]),
        "drop_smallest": list(spec["size_windows"]["drop_smallest"]),
        "diagnostic_large_only": list(spec["size_windows"]["diagnostic_large_only"]),
    }
    tc_by_window: dict[str, dict[str, Any]] = {}
    supports: dict[str, dict[str, Any]] = {}
    nu_bounds: dict[str, dict[str, tuple[float, float]]] = {}
    rows: list[dict[str, Any]] = []
    ratio_rows: list[dict[str, Any]] = []

    for wname, sizes in windows.items():
        wc = _window_curves(curves_all, sizes)
        tcinfo = estimate_joint_tc(wc, branch_center=branch_center, spec=spec, sizes=sizes)
        tc_by_window[wname] = tcinfo
        tc = float(tcinfo["joint_tc"])
        if not np.isfinite(tc):
            continue
        supports[wname] = {}
        nu_bounds[wname] = {}
        for ch in PRIMARY:
            support = build_locked_support(
                wc[ch], tc=tc, nu=1.0, channel=ch, q=0.0,
                x_window=tuple(spec["support"]["x_window"]),
                edge_points=int(spec["support"]["edge_points"]),
                minimum_points_per_size=int(spec["support"]["minimum_points_per_size"]),
            )
            bounds = symmetric_feasible_nu_bounds(
                wc[ch], tc=tc, channel=ch, support=support, spec=spec,
                interpolation=spec["nu_fit"]["primary_interpolation"],
            )
            fit = fit_nu_fixed_tc(
                wc[ch], tc=tc, channel=ch, support=support, spec=spec,
                nu_bounds=bounds, interpolation=spec["nu_fit"]["primary_interpolation"],
            )
            supports[wname][ch] = support
            nu_bounds[wname][ch] = bounds
            row = fit.to_dict()
            row.update({
                "case_label": case.label,
                "p": case.p,
                "window": wname,
                "sizes": ",".join(str(x) for x in sizes),
                "tc_joint": tc,
                "tc_binder": float(tcinfo["binder_roa"]["estimate"]),
                "tc_xi": float(tcinfo["xi_over_L"]["estimate"]),
                "tc_channel_spread": float(tcinfo["channel_spread"]),
                "support_points_by_L": ";".join(f"{L}:{len(support.target_indices[L])}" for L in support.sizes),
            })
            rows.append(row)

        for ch, kind in (("abs_m", "magnetization"), ("chi_abs", "susceptibility")):
            rr = exponent_ratio_at_tc(wc[ch], tc, kind=kind)
            rr.update({
                "case_label": case.label,
                "p": case.p,
                "window": wname,
                "channel": ch,
                "tc_joint": tc,
            })
            ratio_rows.append(rr)

    # Central sensitivity to branch locator half-width.
    tc_sens: list[dict[str, Any]] = []
    full_sizes = windows["full"]
    full_curves = _window_curves(curves_all, full_sizes)
    for half in spec["tc_estimator"]["branch_half_width_sensitivity"]:
        info = estimate_joint_tc(full_curves, branch_center=branch_center, spec=spec, sizes=full_sizes, branch_half_width=float(half))
        tc_sens.append({
            "case_label": case.label, "p": case.p, "branch_half_width": float(half),
            "tc_joint": float(info["joint_tc"]),
            "tc_binder": float(info["binder_roa"]["estimate"]),
            "tc_xi": float(info["xi_over_L"]["estimate"]),
            "tc_channel_spread": float(info["channel_spread"]),
        })

    # Central interpolation and x-window sensitivity with the independently estimated full Tc.
    sens_rows: list[dict[str, Any]] = []
    if "full" in supports and np.isfinite(tc_by_window["full"]["joint_tc"]):
        tc = float(tc_by_window["full"]["joint_tc"])
        for ch in PRIMARY:
            base_support = supports["full"][ch]
            bounds = nu_bounds["full"][ch]
            for interp in [spec["nu_fit"]["primary_interpolation"], *spec["nu_fit"]["sensitivity_interpolations"]]:
                fr = fit_nu_fixed_tc(full_curves[ch], tc=tc, channel=ch, support=base_support, spec=spec, nu_bounds=bounds, interpolation=str(interp))
                sens_rows.append({"case_label": case.label, "p": case.p, "channel": ch, "kind": "interpolation", "setting": str(interp), "nu": fr.nu, "pb": fr.pb, "boundary_hit": fr.boundary_hit})
            for xw in spec["support"]["sensitivity_x_windows"]:
                try:
                    sp = build_locked_support(
                        full_curves[ch], tc=tc, nu=1.0, channel=ch, q=0.0,
                        x_window=tuple(float(x) for x in xw),
                        edge_points=int(spec["support"]["edge_points"]),
                        minimum_points_per_size=int(spec["support"]["minimum_points_per_size"]),
                    )
                    bd = symmetric_feasible_nu_bounds(full_curves[ch], tc=tc, channel=ch, support=sp, spec=spec)
                    fr = fit_nu_fixed_tc(full_curves[ch], tc=tc, channel=ch, support=sp, spec=spec, nu_bounds=bd)
                    sens_rows.append({"case_label": case.label, "p": case.p, "channel": ch, "kind": "x_window", "setting": str(tuple(xw)), "nu": fr.nu, "pb": fr.pb, "boundary_hit": fr.boundary_hit})
                except Exception as exc:
                    sens_rows.append({"case_label": case.label, "p": case.p, "channel": ch, "kind": "x_window", "setting": str(tuple(xw)), "nu": float("nan"), "pb": float("nan"), "boundary_hit": True, "error": f"{type(exc).__name__}: {exc}"})

    return {
        "case_label": case.label,
        "p": case.p,
        "branch_center": branch_center,
        "central_curves": curves_all,
        "tc_by_window": tc_by_window,
        "supports": supports,
        "nu_bounds": nu_bounds,
        "fit_rows": rows,
        "ratio_rows": ratio_rows,
        "tc_sensitivity_rows": tc_sens,
        "sensitivity_rows": sens_rows,
    }


def central_all(cases: dict[str, Any], spec: dict[str, Any]) -> dict[str, Any]:
    built = {label: build_central_case(case, spec) for label, case in cases.items()}
    support_rows: list[dict[str, Any]] = []
    for label, cc in built.items():
        for window, by_ch in cc["supports"].items():
            for ch, sp in by_ch.items():
                counts = {int(L): int(len(sp.target_indices[L])) for L in sp.sizes}
                vals = np.asarray(list(counts.values()), float)
                support_rows.append({
                    "case_label": label, "p": float(cc["p"]), "window": window, "channel": ch,
                    "n_sizes": int(len(counts)), "min_target_points_per_size": int(np.min(vals)),
                    "max_target_points_per_size": int(np.max(vals)),
                    "max_to_min_target_point_ratio": float(np.max(vals) / np.min(vals)),
                    "total_target_points": int(np.sum(vals)),
                    "n_ordered_residuals": int(sp.n_ordered_residuals),
                    "target_points_by_L": ";".join(f"{L}:{counts[L]}" for L in sorted(counts)),
                    "interpretation": "diagnostic only; Bhattacharjee-Seno primary residue remains the published all-point definition and is not posthoc reweighted",
                })
    return {
        "cases": built,
        "fits": pd.DataFrame([r for x in built.values() for r in x["fit_rows"]]),
        "ratios": pd.DataFrame([r for x in built.values() for r in x["ratio_rows"]]),
        "tc_sensitivity": pd.DataFrame([r for x in built.values() for r in x["tc_sensitivity_rows"]]),
        "sensitivity": pd.DataFrame([r for x in built.values() for r in x["sensitivity_rows"]]),
        "support_balance": pd.DataFrame(support_rows),
    }
