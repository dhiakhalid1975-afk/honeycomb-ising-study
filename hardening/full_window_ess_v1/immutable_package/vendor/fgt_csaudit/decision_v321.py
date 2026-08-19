from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

PRIMARY = ("binder_roa", "xi_over_L")


def _row(df: pd.DataFrame, label: str, channel: str | None = None) -> pd.Series:
    q = df.loc[df.case_label == label]
    if channel is not None:
        q = q.loc[q.channel == channel]
    if len(q) != 1:
        raise RuntimeError(f"expected one row for {label}/{channel}, got {len(q)}")
    return q.iloc[0]


def pristine_calibration(summary: dict[str, pd.DataFrame], spec: dict[str, Any], core_audit: dict[str, Any]) -> dict[str, Any]:
    label = str(spec["calibration"]["pristine_case"])
    tc = _row(summary["tc"], label)
    exact = float(spec["exact_honeycomb_ising_tc"])
    tc_ok = bool(float(tc.tc_robustness_low) <= exact <= float(tc.tc_robustness_high))
    nu_checks: dict[str, Any] = {}
    for ch in PRIMARY:
        r = _row(summary["nu"], label, ch)
        ok = bool(float(r.nu_robustness_low) <= 1.0 <= float(r.nu_robustness_high))
        nu_checks[ch] = {
            "pass": ok,
            "bootstrap_ci": [float(r.nu_bootstrap_ci_low), float(r.nu_bootstrap_ci_high)],
            "robustness_envelope": [float(r.nu_robustness_low), float(r.nu_robustness_high)],
            "paired_drift_significant_toward_one": bool(r.paired_drift_significant_toward_one),
        }
    ratio_checks: dict[str, Any] = {}
    for ch in ("abs_m", "chi_abs"):
        r = _row(summary["ratios"], label, ch)
        ref = float(r.ising_reference)
        ok = bool(float(r.ratio_robustness_low) <= ref <= float(r.ratio_robustness_high))
        ratio_checks[ch] = {
            "pass": ok,
            "reference": ref,
            "bootstrap_ci": [float(r.ratio_ci_low), float(r.ratio_ci_high)],
            "robustness_envelope": [float(r.ratio_robustness_low), float(r.ratio_robustness_high)],
            "paired_Lmin_drift_median": float(r.paired_Lmin_drift_median),
        }
    core_ok = bool(core_audit.get("passed", False))
    require_ratio = bool(spec["calibration"].get("require_exponent_ratio_references_inside_robustness_envelopes", True))
    ratio_ok = all(x["pass"] for x in ratio_checks.values()) if require_ratio else True
    passed = bool(tc_ok and all(x["pass"] for x in nu_checks.values()) and ratio_ok and core_ok)
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "exact_tc": exact,
        "tc_pass": tc_ok,
        "tc_bootstrap_ci": [float(tc.tc_bootstrap_ci_low), float(tc.tc_bootstrap_ci_high)],
        "tc_robustness_envelope": [float(tc.tc_robustness_low), float(tc.tc_robustness_high)],
        "nu_checks": nu_checks,
        "exponent_ratio_checks": ratio_checks,
        "exponent_ratio_gate_required": require_ratio,
        "monte_carlo_core_audit_pass": core_ok,
        "note": "Calibration uses exact p=1 theory only as an external validation gate, never to tune diluted-case Tc, nu, beta/nu, or gamma/nu.",
    }


def decide_cases(summary: dict[str, pd.DataFrame], spec: dict[str, Any], calibration: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    labels = list(summary["tc"].case_label)
    max_boundary = float(spec["decision"].get("max_primary_boundary_fraction", 0.10))
    min_valid = float(spec["decision"].get("min_primary_valid_fraction", 0.95))
    for label in labels:
        tc = _row(summary["tc"], label)
        p = float(tc.p)
        prim = [_row(summary["nu"], label, ch) for ch in PRIMARY]
        robust_contains = [float(r.nu_robustness_low) <= 1.0 <= float(r.nu_robustness_high) for r in prim]
        nominal_contains = [float(r.nu_bootstrap_ci_low) <= 1.0 <= float(r.nu_bootstrap_ci_high) for r in prim]
        sides = []
        for r in prim:
            if float(r.nu_robustness_low) > 1.0:
                sides.append("above")
            elif float(r.nu_robustness_high) < 1.0:
                sides.append("below")
            else:
                sides.append("contains")
        drift_toward = [bool(r.paired_drift_significant_toward_one) for r in prim]
        boundary_bad = any(float(r.bootstrap_boundary_fraction) > max_boundary for r in prim)
        validity_bad = any(float(r.bootstrap_valid_fraction) < min_valid for r in prim)
        tc_consistent = bool(tc.tc_channel_ci_overlap)

        reason = ""
        if label == spec["calibration"]["pristine_case"]:
            if calibration["passed"]:
                status = "NU1_COMPATIBLE"
                reason = "exact_pristine_calibration_pass"
            else:
                status = "INCONCLUSIVE_LIMITED_RANGE"
                reason = "exact_pristine_calibration_not_passed"
        elif any(drift_toward):
            status = "INCONCLUSIVE_CORRECTION_DOMINATED"
            reason = "paired_Lmin_drift_significantly_moves_toward_nu1"
        elif all(robust_contains):
            status = "NU1_COMPATIBLE"
            reason = "both_primary_robustness_envelopes_include_nu1"
        elif boundary_bad:
            status = "INCONCLUSIVE_LIMITED_RANGE"
            reason = f"primary_boundary_fraction_exceeds_{max_boundary:.3f}"
        elif validity_bad:
            status = "INCONCLUSIVE_LIMITED_RANGE"
            reason = f"primary_valid_fraction_below_{min_valid:.3f}"
        elif not tc_consistent:
            status = "INCONCLUSIVE_LIMITED_RANGE"
            reason = "binder_and_xi_tc_bootstrap_intervals_do_not_overlap"
        elif not calibration["passed"]:
            status = "INCONCLUSIVE_LIMITED_RANGE"
            reason = "pristine_calibration_veto"
        elif sides[0] == sides[1] and sides[0] in {"above", "below"}:
            status = "EVIDENCE_AGAINST_NU1"
            reason = "both_primary_robustness_envelopes_exclude_nu1_same_side_after_all_vetoes"
        else:
            status = "INCONCLUSIVE_LIMITED_RANGE"
            reason = "primary_channels_not_jointly_decisive"

        rows.append({
            "case_label": label, "p": p, "decision": status, "reason": reason,
            "pristine_calibration_pass": bool(calibration["passed"]),
            "tc_channel_ci_overlap": tc_consistent,
            "binder_robustness_contains_1": robust_contains[0],
            "xi_robustness_contains_1": robust_contains[1],
            "binder_nominal_bootstrap_ci_contains_1": nominal_contains[0],
            "xi_nominal_bootstrap_ci_contains_1": nominal_contains[1],
            "binder_significant_drift_toward_1": drift_toward[0],
            "xi_significant_drift_toward_1": drift_toward[1],
            "primary_boundary_gate_pass": not boundary_bad,
            "primary_validity_gate_pass": not validity_bad,
            "scope": "No decision is proof of a new universality class; logarithmic/finite-size corrections remain a competing explanation.",
        })
    return pd.DataFrame(rows)


def claim_scope_table(decisions: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([
        {"claim": "RG-crossing Tc(p) estimates with quenched-realization bootstrap", "allowed": True, "status": "ALLOWED_WITH_METHOD_LABEL", "note": "Report Binder and xi/L channel estimates and joint estimator; distinguish bootstrap CI from robustness envelope."},
        {"claim": "nu=1 compatibility within current finite-size range", "allowed": True, "status": "ALLOWED_IF_DECISION_SUPPORTS", "note": "Use NU1_COMPATIBLE or INCONCLUSIVE language exactly as generated."},
        {"claim": "finite-size/logarithmic corrections dominate apparent exponent drift", "allowed": True, "status": "ALLOWED_WHEN_PAIRED_DRIFT_GATE_TRIGGERS", "note": "Describe as correction-dominated finite-size evidence, not direct measurement of a logarithmic exponent."},
        {"claim": "nu differs asymptotically from 1", "allowed": False, "status": "PROHIBITED_UNLESS_EVIDENCE_AGAINST_NU1_AND_EXTERNAL_REVIEW", "note": "Even EVIDENCE_AGAINST_NU1 is framework-specific and does not establish a new universality class."},
        {"claim": "new universality class", "allowed": False, "status": "PROHIBITED", "note": "Five sizes over L=40..120 cannot establish a new universality class in the marginal-disorder problem."},
        {"claim": "DFT-parameterized or Kelvin material prediction", "allowed": False, "status": "PROHIBITED", "note": "No material exchange/anisotropy values are injected into the accepted Ising Monte Carlo data."},
    ])
