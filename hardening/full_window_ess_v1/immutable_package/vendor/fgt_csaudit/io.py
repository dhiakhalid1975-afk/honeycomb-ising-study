from __future__ import annotations

from dataclasses import dataclass
import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .provenance import sha256_file




def _nanmean_keep_nan(a: np.ndarray, axis: int = 0) -> np.ndarray:
    """NaN-aware mean without emitting RuntimeWarning when a bootstrap draw is all-NaN.

    All-NaN positions remain NaN and are later treated as failed locked-support candidates.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)

@dataclass(frozen=True)
class CaseData:
    label: str
    p: float
    sizes: tuple[int, ...]
    temperatures: dict[int, np.ndarray]
    values: dict[str, dict[int, np.ndarray]]
    n_realizations: dict[int, int]
    stage: str

    def central_curves(self) -> dict[str, dict[int, tuple[np.ndarray, np.ndarray]]]:
        out: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]] = {k: {} for k in ("abs_m", "chi_abs", "binder_roa", "xi_over_L")}
        for L in self.sizes:
            t = self.temperatures[L]
            out["abs_m"][L] = (t, np.mean(self.values["abs_m"][L], axis=0))
            out["chi_abs"][L] = (t, np.mean(self.values["chi_abs"][L], axis=0))
            m2 = np.mean(self.values["m2"][L], axis=0)
            m4 = np.mean(self.values["m4"][L], axis=0)
            binder = 1.0 - m4 / np.maximum(3.0 * m2 * m2, 1e-300)
            out["binder_roa"][L] = (t, binder)
            xi = self.values["xi_over_L"][L]
            out["xi_over_L"][L] = (t, _nanmean_keep_nan(xi, axis=0))
        return out

    def bootstrap_curves(self, bootstrap_index: int, base_seed: int) -> dict[str, dict[int, tuple[np.ndarray, np.ndarray]]]:
        out: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]] = {k: {} for k in ("abs_m", "chi_abs", "binder_roa", "xi_over_L")}
        # Independent resampling within each L matches the accepted N60 quenched-bootstrap implementation.
        for iL, L in enumerate(self.sizes):
            seed = np.random.SeedSequence([int(base_seed), int(round(self.p * 1000)), int(bootstrap_index), int(iL)])
            rng = np.random.default_rng(seed)
            n = self.n_realizations[L]
            draw = rng.integers(0, n, size=n, endpoint=False)
            t = self.temperatures[L]
            out["abs_m"][L] = (t, np.mean(self.values["abs_m"][L][draw], axis=0))
            out["chi_abs"][L] = (t, np.mean(self.values["chi_abs"][L][draw], axis=0))
            m2 = np.mean(self.values["m2"][L][draw], axis=0)
            m4 = np.mean(self.values["m4"][L][draw], axis=0)
            out["binder_roa"][L] = (t, 1.0 - m4 / np.maximum(3.0 * m2 * m2, 1e-300))
            xi = self.values["xi_over_L"][L][draw]
            out["xi_over_L"][L] = (t, _nanmean_keep_nan(xi, axis=0))
        return out


def _validate_table_shape(df: pd.DataFrame, label: str) -> None:
    required = {"L", "temperature", "realization", "case_label", "p_target", "abs_m", "chi_abs", "m2", "m4", "xi_over_L"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{label}: missing required columns {sorted(missing)}")
    # xi/L can be mathematically undefined for a small number of realization/temperature
    # rows when the finite-size structure-factor estimator is not positive/finite.  The
    # accepted N60 crossing implementation uses finite-value averaging and records failed
    # resamples explicitly.  All other inference columns remain strictly finite.
    strict_cols = sorted((required - {"case_label", "xi_over_L"}))
    if df[strict_cols].isna().any().any():
        bad = [c for c in strict_cols if df[c].isna().any()]
        raise ValueError(f"{label}: NaN found in strict required columns: {bad}")
    xi_counts = df.groupby(["case_label", "L", "temperature"])["xi_over_L"].count()
    if (xi_counts < 2).any():
        bad = xi_counts[xi_counts < 2].head(10)
        raise ValueError(f"{label}: fewer than two finite xi/L realizations at one or more (case,L,T): {bad.to_dict()}")
    dup = df.duplicated(["case_label", "L", "temperature", "realization"])
    if bool(dup.any()):
        raise ValueError(f"{label}: duplicate (case,L,T,realization) rows detected")


def _pivot_complete(group: pd.DataFrame, value: str, L: int) -> tuple[np.ndarray, np.ndarray]:
    piv = group.pivot(index="realization", columns="temperature", values=value).sort_index(axis=0).sort_index(axis=1)
    if value != "xi_over_L" and piv.isna().any().any():
        raise ValueError(f"L={L}, {value}: incomplete realization curves")
    # For xi/L, NaN values are retained at their original realization/temperature
    # coordinates. They are never imputed. Central and bootstrap curves use finite-value
    # means; a bootstrap draw with a non-finite locked target is a failed fit, not a silently
    # dropped residual.
    return piv.columns.to_numpy(float), piv.to_numpy(float)


def build_case(df: pd.DataFrame, label: str, expected_p: float, stage: str) -> CaseData:
    g = df.loc[df["case_label"] == label].copy()
    if g.empty:
        raise ValueError(f"case not found: {label}")
    pvals = g["p_target"].to_numpy(float)
    if not np.allclose(pvals, expected_p, atol=1e-12, rtol=0):
        raise ValueError(f"{label}: p_target does not match locked p={expected_p}")
    sizes = tuple(sorted(int(v) for v in g["L"].unique()))
    if sizes != (40, 60, 80, 100, 120):
        raise ValueError(f"{label}: expected sizes (40,60,80,100,120), got {sizes}")
    temps: dict[int, np.ndarray] = {}
    vals: dict[str, dict[int, np.ndarray]] = {k: {} for k in ("abs_m", "chi_abs", "m2", "m4", "xi_over_L")}
    nr: dict[int, int] = {}
    for L in sizes:
        h = g.loc[g["L"] == L].copy()
        base_t = None
        for col in vals:
            t, a = _pivot_complete(h, col, L)
            if base_t is None:
                base_t = t
            elif not np.array_equal(base_t, t):
                raise ValueError(f"{label}, L={L}: temperature grid mismatch across observables")
            vals[col][L] = a
        assert base_t is not None
        temps[L] = base_t
        nr[L] = int(vals["abs_m"][L].shape[0])
        if nr[L] < 2:
            raise ValueError(f"{label}, L={L}: fewer than two complete realizations")
    return CaseData(label, float(expected_p), sizes, temps, vals, nr, stage)


def resolve_inputs(project_root: Path) -> dict[str, Path]:
    base = project_root / "results" / "publication_strict_phase3"
    paths = {
        "fine_realizations": base / "tables" / "fine_realizations_n60.csv",
        "reference_realizations": base / "tables" / "reference_realizations.csv",
        "tc_table": base / "publication_export" / "TABLE_TC_VS_P.csv",
        "adaptive_decision": base / "manifests" / "adaptive_final_decision.json",
        "quality_gate": base / "quality_gate.json",
    }
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("required N60 inputs are missing:\n" + "\n".join(missing))
    return paths


def validate_project_inputs(project_root: Path, spec: dict[str, Any]) -> tuple[dict[str, Path], dict[str, str]]:
    paths = resolve_inputs(project_root)
    adaptive = json.loads(paths["adaptive_decision"].read_text(encoding="utf-8"))
    quality = json.loads(paths["quality_gate"].read_text(encoding="utf-8"))
    if adaptive.get("status") != spec["required_project_gates"]["adaptive_final_status"]:
        raise RuntimeError(f"adaptive gate is not {spec['required_project_gates']['adaptive_final_status']}")
    if quality.get("status") != spec["required_project_gates"]["quality_gate_status"]:
        raise RuntimeError(f"quality gate is not {spec['required_project_gates']['quality_gate_status']}")

    tc = pd.read_csv(paths["tc_table"])
    for label, ref in spec["legacy_tc_provenance"].items():
        row = tc.loc[tc["case_label"] == label]
        if len(row) != 1:
            raise RuntimeError(f"accepted Tc table does not contain exactly one row for {label}")
        r = row.iloc[0]
        for key, col in (("p", "p"), ("tc", "tc"), ("ci_low", "ci_low"), ("ci_high", "ci_high")):
            if abs(float(r[col]) - float(ref[key])) > 1e-12:
                raise RuntimeError(f"legacy accepted N60 provenance value changed for {label}:{col}; audit refuses to continue")
    exact = float(spec["exact_honeycomb_ising_tc"])
    p1 = spec["legacy_tc_provenance"]["pristine_p100"]
    if not (float(p1["ci_low"]) <= exact <= float(p1["ci_high"])):
        raise RuntimeError("locked pristine accepted CI no longer contains exact honeycomb Ising Tc")

    fine = pd.read_csv(paths["fine_realizations"])
    refdf = pd.read_csv(paths["reference_realizations"])
    _validate_table_shape(fine, "fine_realizations_n60")
    _validate_table_shape(refdf, "reference_realizations")
    hashes = {k: sha256_file(v) for k, v in paths.items()}
    return paths, hashes


def load_all_cases(paths: dict[str, Path], spec: dict[str, Any]) -> dict[str, CaseData]:
    fine = pd.read_csv(paths["fine_realizations"])
    refdf = pd.read_csv(paths["reference_realizations"])
    cases = {
        "random_p080": build_case(fine, "random_p080", 0.80, "fine"),
        "random_p085": build_case(fine, "random_p085", 0.85, "fine"),
        "random_p090": build_case(fine, "random_p090", 0.90, "fine"),
        "pristine_p100": build_case(refdf, "pristine_p100", 1.00, "reference"),
    }
    return cases


def input_data_quality_table(paths: dict[str, Path]) -> pd.DataFrame:
    """Summarize finite-value coverage without modifying or imputing source data."""
    rows: list[dict[str, Any]] = []
    for source_key in ("fine_realizations", "reference_realizations"):
        df = pd.read_csv(paths[source_key])
        for (label, L), g in df.groupby(["case_label", "L"], sort=True):
            xi = g["xi_over_L"]
            finite_per_t = g.groupby("temperature")["xi_over_L"].count()
            rows.append({
                "source": source_key,
                "case_label": str(label),
                "L": int(L),
                "n_rows": int(len(g)),
                "n_temperatures": int(g["temperature"].nunique()),
                "n_realizations": int(g["realization"].nunique()),
                "xi_over_L_nonfinite_rows": int((~np.isfinite(xi.to_numpy(float))).sum()),
                "xi_over_L_nonfinite_fraction": float((~np.isfinite(xi.to_numpy(float))).mean()),
                "min_finite_xi_realizations_per_temperature": int(finite_per_t.min()),
                "strict_observable_nonfinite_rows": int(sum((~np.isfinite(g[c].to_numpy(float))).sum() for c in ("abs_m", "chi_abs", "m2", "m4"))),
                "handling": "xi/L nonfinite values retained; finite-value mean only; locked-target nonfinite bootstrap fits fail closed; no imputation",
            })
    return pd.DataFrame(rows)
