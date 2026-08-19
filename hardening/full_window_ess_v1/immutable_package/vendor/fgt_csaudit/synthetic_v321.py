from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .correction_fit import exponent_ratio_at_tc, fit_nu_fixed_tc, restrict_curves, symmetric_feasible_nu_bounds
from .pb import build_locked_support
from .provenance import atomic_write_json
from .rg_tc import estimate_joint_tc

SIZES = (40, 60, 80, 100, 120)


def _synthetic_curves(*, tc: float, nu: float, correction_amp: float = 0.0, omega: float = 1.0) -> dict[str, dict[int, tuple[np.ndarray, np.ndarray]]]:
    out = {k: {} for k in ("binder_roa", "xi_over_L", "abs_m", "chi_abs")}
    for L in SIZES:
        T = np.linspace(tc - 0.075, tc + 0.075, 31)
        x = (T - tc) / tc * L ** (1.0 / nu)
        corr = correction_amp * L ** (-omega)
        binder = 0.610 - 0.105 * np.tanh(x / 2.2) + corr
        xi = 0.920 - 0.170 * np.tanh(x / 2.0) + 1.6 * corr
        fm = 0.80 + 0.10 * np.tanh(-x / 2.0)
        fc = 0.95 + 0.45 * np.exp(-0.25 * x * x)
        absm = L ** (-0.125) * fm
        chi = L ** 1.75 * fc
        out["binder_roa"][L] = (T, binder)
        out["xi_over_L"][L] = (T, xi)
        out["abs_m"][L] = (T, absm)
        out["chi_abs"][L] = (T, chi)
    return out


def _fit(curves: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]], spec: dict[str, Any], branch: float, sizes: tuple[int, ...]) -> dict[str, float]:
    tcinfo = estimate_joint_tc(curves, branch_center=branch, spec=spec, sizes=sizes)
    tc = float(tcinfo["joint_tc"])
    res: dict[str, float] = {"tc": tc}
    for ch in ("binder_roa", "xi_over_L"):
        c = restrict_curves(curves[ch], sizes)
        sp = build_locked_support(c, tc=tc, nu=1.0, channel=ch, q=0.0, x_window=tuple(spec["support"]["x_window"]), edge_points=int(spec["support"]["edge_points"]), minimum_points_per_size=int(spec["support"]["minimum_points_per_size"]))
        bd = symmetric_feasible_nu_bounds(c, tc=tc, channel=ch, support=sp, spec=spec)
        fr = fit_nu_fixed_tc(c, tc=tc, channel=ch, support=sp, spec=spec, nu_bounds=bd)
        res[f"nu_{ch}"] = float(fr.nu)
    return res


def run_synthetic_challenge(out_dir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def add(name: str, ok: bool, detail: Any) -> None:
        checks.append({"check": name, "pass": bool(ok), "detail": detail})

    # Ideal recovery.
    tc0 = 1.2345
    ideal = _synthetic_curves(tc=tc0, nu=1.0)
    f = _fit(ideal, spec, tc0, SIZES)
    add("ideal_tc_recovery", abs(f["tc"] - tc0) < 0.002, f)
    add("ideal_binder_nu_recovery", abs(f["nu_binder_roa"] - 1.0) < 0.02, f["nu_binder_roa"])
    add("ideal_xi_nu_recovery", abs(f["nu_xi_over_L"] - 1.0) < 0.02, f["nu_xi_over_L"])
    rbeta = exponent_ratio_at_tc(ideal["abs_m"], f["tc"], kind="magnetization")
    rgamma = exponent_ratio_at_tc(ideal["chi_abs"], f["tc"], kind="susceptibility")
    add("ideal_beta_over_nu_recovery", bool(rbeta.get("success")) and abs(float(rbeta["ratio"]) - 0.125) < 0.01, rbeta)
    add("ideal_gamma_over_nu_recovery", bool(rgamma.get("success")) and abs(float(rgamma["ratio"]) - 1.75) < 0.03, rgamma)

    # Correction challenge: current-size effective exponent should move toward 1 after dropping L=40.
    corr = _synthetic_curves(tc=tc0, nu=1.0, correction_amp=0.65, omega=0.7)
    ff = _fit(corr, spec, tc0, SIZES)
    fd = _fit(corr, spec, tc0, SIZES[1:])
    add("correction_challenge_binder_drift_toward_one", abs(fd["nu_binder_roa"] - 1.0) <= abs(ff["nu_binder_roa"] - 1.0) + 1e-6, {"full": ff, "drop": fd})
    add("correction_challenge_xi_drift_toward_one", abs(fd["nu_xi_over_L"] - 1.0) <= abs(ff["nu_xi_over_L"] - 1.0) + 1e-6, {"full": ff, "drop": fd})

    # True shifted exponent without corrections should be comparatively stable with Lmin.
    shifted = _synthetic_curves(tc=tc0, nu=1.12)
    sf = _fit(shifted, spec, tc0, SIZES)
    sd = _fit(shifted, spec, tc0, SIZES[1:])
    add("true_shifted_nu_binder_recovered", abs(sf["nu_binder_roa"] - 1.12) < 0.04, {"full": sf, "drop": sd})
    add("true_shifted_nu_xi_recovered", abs(sf["nu_xi_over_L"] - 1.12) < 0.04, {"full": sf, "drop": sd})
    add("true_shifted_nu_Lmin_stability", max(abs(sd["nu_binder_roa"] - sf["nu_binder_roa"]), abs(sd["nu_xi_over_L"] - sf["nu_xi_over_L"])) < 0.04, {"full": sf, "drop": sd})

    passed = all(x["pass"] for x in checks)
    result = {"status": "PASS" if passed else "FAIL", "passed": passed, "n_checks": len(checks), "n_failed": sum(not x["pass"] for x in checks), "checks": checks}
    d = Path(out_dir) / "synthetic"
    d.mkdir(parents=True, exist_ok=True)
    atomic_write_json(d / "SYNTHETIC_CHALLENGE_RESULT.json", result)
    if not passed:
        raise RuntimeError("synthetic correction-aware challenge failed")
    return result
