from __future__ import annotations

from typing import Any
import pandas as pd

PRIMARY = ("binder_roa", "xi_over_L")


def _one(df: pd.DataFrame, label: str, channel: str | None = None) -> pd.Series:
    q = df.loc[df.case_label == label]
    if channel is not None:
        q = q.loc[q.channel == channel]
    if len(q) != 1:
        raise RuntimeError(f"expected one row for {label}/{channel}, got {len(q)}")
    return q.iloc[0]


def _as_bool(v: Any) -> bool:
    if isinstance(v, str):
        return v.strip().lower() in {"true", "1", "yes"}
    return bool(v)


def corrected_decisions(
    tc: pd.DataFrame,
    nu: pd.DataFrame,
    *,
    pristine_calibration_pass: bool,
    pristine_case: str = "pristine_p100",
    max_primary_boundary_fraction: float = 0.10,
    min_primary_valid_fraction: float = 0.95,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in tc.case_label.tolist():
        tr = _one(tc, label)
        prim = [_one(nu, label, ch) for ch in PRIMARY]
        robust_contains = [float(r.nu_robustness_low) <= 1.0 <= float(r.nu_robustness_high) for r in prim]
        nominal_contains = [float(r.nu_bootstrap_ci_low) <= 1.0 <= float(r.nu_bootstrap_ci_high) for r in prim]
        drift_toward = [_as_bool(r.paired_drift_significant_toward_one) for r in prim]
        boundary_bad = any(float(r.bootstrap_boundary_fraction) > max_primary_boundary_fraction for r in prim)
        validity_bad = any(float(r.bootstrap_valid_fraction) < min_primary_valid_fraction for r in prim)
        tc_consistent = _as_bool(tr.tc_channel_ci_overlap)
        sides: list[str] = []
        for r in prim:
            lo, hi = float(r.nu_robustness_low), float(r.nu_robustness_high)
            sides.append("above" if lo > 1.0 else "below" if hi < 1.0 else "contains")

        if label == pristine_case:
            if pristine_calibration_pass:
                status, reason = "NU1_COMPATIBLE", "exact_pristine_calibration_pass"
            else:
                status, reason = "INCONCLUSIVE_LIMITED_RANGE", "exact_pristine_calibration_not_passed"
        elif boundary_bad:
            status, reason = "INCONCLUSIVE_LIMITED_RANGE", f"primary_boundary_fraction_exceeds_{max_primary_boundary_fraction:.3f}"
        elif validity_bad:
            status, reason = "INCONCLUSIVE_LIMITED_RANGE", f"primary_valid_fraction_below_{min_primary_valid_fraction:.3f}"
        elif not tc_consistent:
            status, reason = "INCONCLUSIVE_LIMITED_RANGE", "binder_and_xi_tc_bootstrap_intervals_do_not_overlap"
        elif not pristine_calibration_pass:
            status, reason = "INCONCLUSIVE_LIMITED_RANGE", "pristine_calibration_veto"
        elif any(drift_toward):
            status, reason = "INCONCLUSIVE_CORRECTION_DOMINATED", "paired_Lmin_drift_significantly_moves_toward_nu1"
        elif all(robust_contains):
            status, reason = "NU1_COMPATIBLE", "both_primary_robustness_envelopes_include_nu1_after_all_vetoes"
        elif sides[0] == sides[1] and sides[0] in {"above", "below"}:
            status, reason = "EVIDENCE_AGAINST_NU1", "both_primary_robustness_envelopes_exclude_nu1_same_side_after_all_vetoes"
        else:
            status, reason = "INCONCLUSIVE_LIMITED_RANGE", "primary_channels_not_jointly_decisive"

        rows.append({
            "case_label": label,
            "p": float(tr.p),
            "decision": status,
            "reason": reason,
            "pristine_calibration_pass": bool(pristine_calibration_pass),
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


def claim_scope_table() -> pd.DataFrame:
    return pd.DataFrame([
        {"claim": "RG-crossing Tc(p) estimates with quenched-realization bootstrap", "allowed": True, "status": "ALLOWED_WITH_METHOD_LABEL", "note": "Report Binder and xi/L channel estimates and joint estimator; distinguish bootstrap CI from robustness envelope."},
        {"claim": "nu=1 compatibility within current finite-size range", "allowed": True, "status": "ALLOWED_IF_DECISION_SUPPORTS", "note": "Use NU1_COMPATIBLE or INCONCLUSIVE language exactly as generated."},
        {"claim": "finite-size/logarithmic corrections dominate apparent exponent drift", "allowed": True, "status": "ALLOWED_WHEN_PAIRED_DRIFT_GATE_TRIGGERS", "note": "Describe as correction-dominated finite-size evidence, not direct measurement of a logarithmic exponent."},
        {"claim": "nu differs asymptotically from 1", "allowed": False, "status": "PROHIBITED_UNLESS_EVIDENCE_AGAINST_NU1_AND_EXTERNAL_REVIEW", "note": "Even EVIDENCE_AGAINST_NU1 is framework-specific and does not establish a new universality class."},
        {"claim": "new universality class", "allowed": False, "status": "PROHIBITED", "note": "Five sizes over L=40..120 cannot establish a new universality class in the marginal-disorder problem."},
        {"claim": "DFT-parameterized or Kelvin material prediction", "allowed": False, "status": "PROHIBITED", "note": "No material exchange/anisotropy values are injected into the accepted Ising Monte Carlo data."},
    ])
