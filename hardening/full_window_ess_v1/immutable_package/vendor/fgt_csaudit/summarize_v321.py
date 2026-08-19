from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .analysis_v321 import PRIMARY


def _q(a: pd.Series | np.ndarray, p: float) -> float:
    x = np.asarray(a, float)
    x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if len(x) else float("nan")


def _median(a: pd.Series | np.ndarray) -> float:
    x = np.asarray(a, float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else float("nan")


def _summary(a: pd.Series | np.ndarray) -> dict[str, float | int]:
    x = np.asarray(a, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {"n": 0, "median": float("nan"), "mean": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    return {
        "n": int(len(x)), "median": float(np.median(x)), "mean": float(np.mean(x)),
        "ci_low": float(np.quantile(x, 0.025)), "ci_high": float(np.quantile(x, 0.975)),
    }


def _central_nu(central_fits: pd.DataFrame, label: str, ch: str, window: str) -> float:
    r = central_fits.loc[(central_fits.case_label == label) & (central_fits.channel == ch) & (central_fits.window == window)]
    return float(r.iloc[0].nu) if len(r) == 1 else float("nan")


def summarize_production(boot: pd.DataFrame, central: dict[str, Any], spec: dict[str, Any]) -> dict[str, pd.DataFrame]:
    labels = list(central["cases"])
    central_fits = central["fits"]
    sens = central["sensitivity"]
    tc_sens = central["tc_sensitivity"]

    tc_rows: list[dict[str, Any]] = []
    nu_rows: list[dict[str, Any]] = []
    drift_rows: list[dict[str, Any]] = []
    ratio_rows: list[dict[str, Any]] = []

    for label in labels:
        g = boot.loc[boot.case_label == label].copy()
        p = float(central["cases"][label]["p"])
        cc = central["cases"][label]
        tc_full_c = float(cc["tc_by_window"]["full"]["joint_tc"])
        tc_drop_c = float(cc["tc_by_window"]["drop_smallest"]["joint_tc"])
        tc_b_c = float(cc["tc_by_window"]["full"]["binder_roa"]["estimate"])
        tc_x_c = float(cc["tc_by_window"]["full"]["xi_over_L"]["estimate"])

        sj = _summary(g.tc_full_joint); sb = _summary(g.tc_full_binder); sx = _summary(g.tc_full_xi)
        branch_vals = tc_sens.loc[tc_sens.case_label == label, "tc_joint"].to_numpy(float)
        branch_delta = float(np.nanmax(np.abs(branch_vals - tc_full_c))) if len(branch_vals) and np.isfinite(branch_vals).any() else 0.0
        tc_sys = max(
            abs(tc_drop_c - tc_full_c) if np.isfinite(tc_drop_c) else 0.0,
            0.5 * abs(tc_b_c - tc_x_c) if np.isfinite(tc_b_c) and np.isfinite(tc_x_c) else 0.0,
            branch_delta,
        )
        tc_low = float(sj["ci_low"] - tc_sys) if np.isfinite(float(sj["ci_low"])) else float("nan")
        tc_high = float(sj["ci_high"] + tc_sys) if np.isfinite(float(sj["ci_high"])) else float("nan")
        tc_overlap = bool(max(float(sb["ci_low"]), float(sx["ci_low"])) <= min(float(sb["ci_high"]), float(sx["ci_high"]))) if all(np.isfinite([sb["ci_low"], sb["ci_high"], sx["ci_low"], sx["ci_high"]])) else False
        tc_rows.append({
            "case_label": label, "p": p, "tc_central_joint": tc_full_c,
            "tc_bootstrap_median": sj["median"], "tc_bootstrap_ci_low": sj["ci_low"], "tc_bootstrap_ci_high": sj["ci_high"],
            "tc_binder_median": sb["median"], "tc_binder_ci_low": sb["ci_low"], "tc_binder_ci_high": sb["ci_high"],
            "tc_xi_median": sx["median"], "tc_xi_ci_low": sx["ci_low"], "tc_xi_ci_high": sx["ci_high"],
            "tc_channel_ci_overlap": tc_overlap, "tc_systematic_radius": tc_sys,
            "tc_robustness_low": tc_low, "tc_robustness_high": tc_high,
            "robustness_note": "bootstrap percentile interval expanded by max central Lmin/channel/branch sensitivity; not a confidence interval",
        })

        for ch in PRIMARY:
            full_col = f"nu_{ch}_full"; drop_col = f"nu_{ch}_drop"
            sfull = _summary(g[full_col]); sdrop = _summary(g[drop_col])
            drift = g[drop_col].to_numpy(float) - g[full_col].to_numpy(float)
            sd = _summary(drift)
            central_full = _central_nu(central_fits, label, ch, "full")
            central_drop = _central_nu(central_fits, label, ch, "drop_smallest")
            # Tc channel-choice systematic from the same bootstrap draws.
            d_b = np.abs(g[f"nu_{ch}_binderTc"].to_numpy(float) - g[full_col].to_numpy(float))
            d_x = np.abs(g[f"nu_{ch}_xiTc"].to_numpy(float) - g[full_col].to_numpy(float))
            tc_choice = max(_median(d_b), _median(d_x))
            # Central interpolation / x-window sensitivity relative to primary central fit.
            ss = sens.loc[(sens.case_label == label) & (sens.channel == ch)].copy()
            deviations = np.abs(ss.nu.to_numpy(float) - central_full)
            deviations = deviations[np.isfinite(deviations)]
            central_sens = float(np.max(deviations)) if len(deviations) else 0.0
            drift_mag = abs(float(sd["median"])) if np.isfinite(float(sd["median"])) else 0.0
            sys_radius = max(drift_mag, tc_choice if np.isfinite(tc_choice) else 0.0, central_sens)
            rob_low = float(sfull["ci_low"] - sys_radius) if np.isfinite(float(sfull["ci_low"])) else float("nan")
            rob_high = float(sfull["ci_high"] + sys_radius) if np.isfinite(float(sfull["ci_high"])) else float("nan")
            bfrac = float(np.mean(g[f"boundary_{ch}_full"].astype(bool)))
            validfrac = float(np.mean(g[f"valid_{ch}_full"].astype(bool)))
            full_med = float(sfull["median"])
            toward = False
            significant_toward = False
            if np.isfinite(full_med) and np.isfinite(float(sd["ci_low"])) and np.isfinite(float(sd["ci_high"])):
                if full_med > 1.0:
                    toward = float(sd["median"]) < 0.0
                    significant_toward = float(sd["ci_high"]) < 0.0
                elif full_med < 1.0:
                    toward = float(sd["median"]) > 0.0
                    significant_toward = float(sd["ci_low"]) > 0.0
            nu_rows.append({
                "case_label": label, "p": p, "channel": ch,
                "nu_central_full": central_full, "nu_central_drop_smallest": central_drop,
                "nu_bootstrap_median": sfull["median"], "nu_bootstrap_ci_low": sfull["ci_low"], "nu_bootstrap_ci_high": sfull["ci_high"],
                "nu_drop_median": sdrop["median"], "nu_drop_ci_low": sdrop["ci_low"], "nu_drop_ci_high": sdrop["ci_high"],
                "paired_drift_median": sd["median"], "paired_drift_ci_low": sd["ci_low"], "paired_drift_ci_high": sd["ci_high"],
                "paired_drift_toward_one": toward, "paired_drift_significant_toward_one": significant_toward,
                "tc_choice_systematic": tc_choice, "central_interpolation_xwindow_systematic": central_sens,
                "systematic_radius": sys_radius, "nu_robustness_low": rob_low, "nu_robustness_high": rob_high,
                "bootstrap_boundary_fraction": bfrac, "bootstrap_valid_fraction": validfrac,
                "robustness_note": "bootstrap percentile CI plus separate conservative systematic envelope; envelope is not a confidence interval",
            })
            drift_rows.append({
                "case_label": label, "p": p, "channel": ch,
                "nu_full_median": sfull["median"], "nu_drop_median": sdrop["median"],
                "delta_nu_drop_minus_full_median": sd["median"], "delta_ci_low": sd["ci_low"], "delta_ci_high": sd["ci_high"],
                "toward_one": toward, "significant_toward_one": significant_toward,
            })

        for ch, ref in (("abs_m", float(spec["reference_exponents"]["beta_over_nu"])), ("chi_abs", float(spec["reference_exponents"]["gamma_over_nu"]))):
            full = g[f"ratio_{ch}_full"].to_numpy(float)
            drop = g[f"ratio_{ch}_drop"].to_numpy(float)
            s1 = _summary(full); s2 = _summary(drop)
            paired = drop - full
            sp = _summary(paired)
            # This is a deliberately conservative finite-size robustness envelope, not a CI.
            # It expands the quenched percentile interval by the median paired Lmin drift.
            ratio_sys = abs(float(sp["median"])) if np.isfinite(float(sp["median"])) else 0.0
            ratio_rob_low = float(s1["ci_low"] - ratio_sys) if np.isfinite(float(s1["ci_low"])) else float("nan")
            ratio_rob_high = float(s1["ci_high"] + ratio_sys) if np.isfinite(float(s1["ci_high"])) else float("nan")
            ratio_rows.append({
                "case_label": label, "p": p, "channel": ch, "ising_reference": ref,
                "ratio_median": s1["median"], "ratio_ci_low": s1["ci_low"], "ratio_ci_high": s1["ci_high"],
                "ratio_drop_median": s2["median"], "ratio_drop_ci_low": s2["ci_low"], "ratio_drop_ci_high": s2["ci_high"],
                "paired_Lmin_drift_median": sp["median"], "paired_Lmin_drift_ci_low": sp["ci_low"], "paired_Lmin_drift_ci_high": sp["ci_high"],
                "ratio_systematic_radius": ratio_sys, "ratio_robustness_low": ratio_rob_low, "ratio_robustness_high": ratio_rob_high,
                "reference_inside_robustness_envelope": bool(np.isfinite(ratio_rob_low) and ratio_rob_low <= ref <= ratio_rob_high),
                "interpretation_scope": "consistency diagnostic at independently estimated RG-crossing Tc; robustness envelope is not a confidence interval and is not a universality proof",
            })

    tc_df = pd.DataFrame(tc_rows)
    nu_df = pd.DataFrame(nu_rows)
    drift_df = pd.DataFrame(drift_rows)
    ratio_df = pd.DataFrame(ratio_rows)

    # Channel-difference diagnostic, paired within each bootstrap draw.
    delta_rows: list[dict[str, Any]] = []
    for label in labels:
        g = boot.loc[boot.case_label == label]
        d = g["nu_binder_roa_full"].to_numpy(float) - g["nu_xi_over_L_full"].to_numpy(float)
        s = _summary(d)
        delta_rows.append({"case_label": label, "p": float(central["cases"][label]["p"]), "delta_nu_binder_minus_xi_median": s["median"], "delta_ci_low": s["ci_low"], "delta_ci_high": s["ci_high"]})
    delta_df = pd.DataFrame(delta_rows)
    return {"tc": tc_df, "nu": nu_df, "drift": drift_df, "ratios": ratio_df, "channel_delta": delta_df}
