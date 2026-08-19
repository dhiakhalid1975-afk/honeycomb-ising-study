from __future__ import annotations
import pandas as pd
from decision_v321_1 import corrected_decisions


def _tables(boundary=(0.20, 0.00), robust=((0.8, 1.2), (0.8, 1.2)), drift=(False, False)):
    tc = pd.DataFrame([{"case_label":"x","p":0.8,"tc_channel_ci_overlap":True}])
    nu = pd.DataFrame([
        {"case_label":"x","channel":"binder_roa","nu_robustness_low":robust[0][0],"nu_robustness_high":robust[0][1],"nu_bootstrap_ci_low":0.9,"nu_bootstrap_ci_high":1.1,"paired_drift_significant_toward_one":drift[0],"bootstrap_boundary_fraction":boundary[0],"bootstrap_valid_fraction":1.0},
        {"case_label":"x","channel":"xi_over_L","nu_robustness_low":robust[1][0],"nu_robustness_high":robust[1][1],"nu_bootstrap_ci_low":0.9,"nu_bootstrap_ci_high":1.1,"paired_drift_significant_toward_one":drift[1],"bootstrap_boundary_fraction":boundary[1],"bootstrap_valid_fraction":1.0},
    ])
    return tc, nu


def main() -> int:
    tc, nu = _tables()
    d = corrected_decisions(tc, nu, pristine_calibration_pass=True).iloc[0]
    assert d.decision == "INCONCLUSIVE_LIMITED_RANGE" and not bool(d.primary_boundary_gate_pass)
    tc, nu = _tables(boundary=(0.0,0.0))
    assert corrected_decisions(tc, nu, pristine_calibration_pass=True).iloc[0].decision == "NU1_COMPATIBLE"
    tc, nu = _tables(boundary=(0.0,0.0), drift=(True,False))
    assert corrected_decisions(tc, nu, pristine_calibration_pass=True).iloc[0].decision == "INCONCLUSIVE_CORRECTION_DOMINATED"
    tc, nu = _tables(boundary=(0.0,0.0), robust=((1.1,1.3),(1.05,1.2)))
    assert corrected_decisions(tc, nu, pristine_calibration_pass=True).iloc[0].decision == "EVIDENCE_AGAINST_NU1"
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
