from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .analysis_v321 import central_all
from .bootstrap_v321 import collect_bootstrap, run_bootstrap_chunks
from .config import load_spec, package_root, spec_hash
from .core_audit import audit_monte_carlo_core
from .decision_v321 import claim_scope_table, decide_cases, pristine_calibration
from .figures_v321 import build_figures
from .io import input_data_quality_table, load_all_cases, validate_project_inputs
from .parallel import backend_name
from .provenance import (
    atomic_write_csv, atomic_write_json, atomic_write_text, environment_record,
    hash_tree, sha256_file, stable_json_hash, utc_now,
)
from .runtime import checkpoint_workspace, checkpoint_workspace_record
from .summarize_v321 import summarize_production
from .synthetic_v321 import run_synthetic_challenge

CASE_LABELS = ["random_p080", "random_p085", "random_p090", "pristine_p100"]


def output_dir(project_root: Path, spec: dict[str, Any]) -> Path:
    return Path(project_root) / Path(spec["write_subdir"])


def _source_hashes() -> dict[str, str]:
    root = package_root()
    out: dict[str, str] = {}
    for sub in ("src", "configs", "literature", "scripts", "tests"):
        d = root / sub
        if d.exists():
            for rel, h in hash_tree(d, suffixes=(".py", ".json", ".md", ".csv", ".toml", ".txt", ".cmd")).items():
                out[f"{sub}/{rel}"] = h
    for name in ("pyproject.toml", "REQUIREMENTS.txt", "DECISIONS_LOCKED.md", "README_AR.md", "README_EN.md", "LICENSE", "VERSION.txt"):
        q = root / name
        if q.exists():
            out[name] = sha256_file(q)
    for q in sorted(root.glob("*.cmd")):
        out[q.name] = sha256_file(q)
    return out


def _read_state(out: Path) -> dict[str, Any]:
    p = out / "RUN_STATE.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _write_state(out: Path, state: dict[str, Any]) -> None:
    atomic_write_json(out / "RUN_STATE.json", state | {"last_update_utc": utc_now()})


def initialize_run(project_root: Path, spec_path: Path | None = None) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    project_root = Path(project_root).resolve()
    spec = load_spec(spec_path)
    paths, input_hashes = validate_project_inputs(project_root, spec)
    core = audit_monte_carlo_core(project_root)
    if not core["passed"]:
        raise RuntimeError("original Monte Carlo core audit failed")
    out = output_dir(project_root, spec)
    out.mkdir(parents=True, exist_ok=True)
    workspace = checkpoint_workspace_record(project_root, out, spec)
    src_hashes = _source_hashes()
    locked = {
        "spec_name": spec["spec_name"], "spec_version": spec["spec_version"],
        "spec_hash": spec_hash(spec), "project_root": str(project_root),
        "input_hashes": input_hashes, "source_hashes": src_hashes,
        "original_core_source_hashes": core["source_hashes"],
        "scientific_scope": spec["scientific_scope"], "executor_backend": backend_name(),
    }
    run_sig = stable_json_hash(locked)
    manifest = out / "manifests" / "RUN_MANIFEST.json"
    if manifest.exists():
        old = json.loads(manifest.read_text(encoding="utf-8"))
        if old.get("run_signature") != run_sig:
            raise RuntimeError("existing v3.2.1 output belongs to a different source/spec/input state; do not mix runs")
    else:
        atomic_write_json(manifest, locked | {"run_signature": run_sig, "created_utc": utc_now()})
    atomic_write_json(out / "manifests" / "MONTE_CARLO_CORE_AUDIT.json", core)
    atomic_write_json(out / "manifests" / "INPUT_MANIFEST.json", {"created_utc": utc_now(), "project_root": str(project_root), "files": {k: {"path": str(paths[k]), "sha256": v} for k, v in input_hashes.items()}})
    atomic_write_json(out / "manifests" / "ENVIRONMENT.json", environment_record())
    atomic_write_json(out / "manifests" / "RUNTIME_WORKSPACE.json", workspace)
    atomic_write_json(out / "manifests" / "SPEC_LOCK_USED.json", spec)
    atomic_write_text(out / "manifests" / "RUN_SIGNATURE.txt", run_sig + "\n")
    atomic_write_csv(out / "tables" / "TABLE_INPUT_DATA_QUALITY.csv", input_data_quality_table(paths))
    st = _read_state(out); st.update({"initialized": True, "run_signature": run_sig}); _write_state(out, st)
    return spec, out, {"paths": paths, "input_hashes": input_hashes, "run_signature": run_sig, "core_audit": core, "runtime_workspace": workspace}


def _central_support_payload(central: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for label, cc in central["cases"].items():
        out[label] = {
            "branch_center": cc["branch_center"],
            "tc_by_window": cc["tc_by_window"],
            "nu_bounds": {w: {ch: list(b) for ch, b in d.items()} for w, d in cc["nu_bounds"].items()},
            "supports": {w: {ch: sp.to_jsonable() for ch, sp in d.items()} for w, d in cc["supports"].items()},
        }
    return out


def ensure_central(project_root: Path, spec_path: Path | None = None) -> dict[str, Any]:
    spec, out, ctx = initialize_run(project_root, spec_path)
    cases = load_all_cases(ctx["paths"], spec)
    print("[central-v3.2.1] independent RG-invariant Tc + symmetric scalar nu fits", flush=True)
    central = central_all(cases, spec)
    tables = out / "tables"; tables.mkdir(parents=True, exist_ok=True)
    atomic_write_csv(tables / "TABLE_CENTRAL_RG_TC_AND_NU.csv", central["fits"])
    atomic_write_csv(tables / "TABLE_CENTRAL_EXPONENT_RATIO_DIAGNOSTICS.csv", central["ratios"])
    atomic_write_csv(tables / "TABLE_TC_BRANCH_SENSITIVITY.csv", central["tc_sensitivity"])
    atomic_write_csv(tables / "TABLE_NU_INTERPOLATION_XWINDOW_SENSITIVITY.csv", central["sensitivity"])
    atomic_write_csv(tables / "TABLE_SUPPORT_BALANCE_AUDIT.csv", central["support_balance"])
    atomic_write_json(tables / "LOCKED_SUPPORT_AND_SYMMETRIC_NU_BOUNDS.json", _central_support_payload(central))
    st = _read_state(out); st["central"] = {"status": "PASS", "table_sha256": sha256_file(tables / "TABLE_CENTRAL_RG_TC_AND_NU.csv")}; _write_state(out, st)
    return {"spec": spec, "out": out, "ctx": ctx, "cases_raw": cases, "central": central}


def run_synthetic(project_root: Path, spec_path: Path | None = None, workers: int = 4) -> dict[str, Any]:
    spec, out, _ = initialize_run(project_root, spec_path)
    _ = workers
    result = run_synthetic_challenge(out, spec)
    st = _read_state(out); st["synthetic"] = result; _write_state(out, st)
    return result


def _require_synthetic(out: Path) -> None:
    p = out / "synthetic" / "SYNTHETIC_CHALLENGE_RESULT.json"
    if not p.exists() or not json.loads(p.read_text(encoding="utf-8")).get("passed", False):
        raise RuntimeError("v3.2.1 synthetic challenge is not PASS")


def run_bootstrap_convergence(project_root: Path, spec_path: Path | None = None, workers: int = 4, force: bool = False) -> dict[str, Any]:
    c = ensure_central(project_root, spec_path)
    spec, out = c["spec"], c["out"]
    _require_synthetic(out)
    n_small = int(spec["bootstrap"]["n_small"]); n_large = int(spec["bootstrap"]["n_large"])
    namespace = f"production_v321_n{n_large}"
    info = run_bootstrap_chunks(project_root=project_root, out_dir=out, spec_path=spec_path, spec=spec, central=c["central"], n_bootstrap=n_large, workers=workers, cases=["random_p080"], namespace=namespace, force=force)
    df = collect_bootstrap(out, checkpoint_root=Path(info["checkpoint_root"]), namespace=namespace, labels=["random_p080"], n_bootstrap=n_large)
    metrics = ["tc_full_joint", "nu_binder_roa_full", "nu_xi_over_L_full"]
    checks: list[dict[str, Any]] = []
    tolfrac = float(spec["bootstrap"]["ci_endpoint_tolerance_fraction_of_large_ci_width"])
    for col in metrics:
        a = df.loc[df.bootstrap_index < n_small, col].to_numpy(float); a = a[np.isfinite(a)]
        b = df.loc[df.bootstrap_index < n_large, col].to_numpy(float); b = b[np.isfinite(b)]
        qa = np.quantile(a, [0.025, 0.975]); qb = np.quantile(b, [0.025, 0.975])
        width = float(qb[1] - qb[0]); tol = tolfrac * width
        shift_lo = abs(float(qa[0] - qb[0])); shift_hi = abs(float(qa[1] - qb[1]))
        checks.append({"metric": col, "small_ci": qa.tolist(), "large_ci": qb.tolist(), "large_width": width, "endpoint_tolerance": tol, "shift_low": shift_lo, "shift_high": shift_hi, "pass": bool(shift_lo <= tol and shift_hi <= tol)})
    # Paired drift convergence is also mandatory.
    for ch in ("binder_roa", "xi_over_L"):
        col = f"drift_{ch}"
        work = df.copy(); work[col] = work[f"nu_{ch}_drop"] - work[f"nu_{ch}_full"]
        a = work.loc[work.bootstrap_index < n_small, col].to_numpy(float); a = a[np.isfinite(a)]
        b = work.loc[work.bootstrap_index < n_large, col].to_numpy(float); b = b[np.isfinite(b)]
        qa = np.quantile(a, [0.025, 0.975]); qb = np.quantile(b, [0.025, 0.975])
        width = float(qb[1] - qb[0]); tol = tolfrac * width
        sl = abs(float(qa[0] - qb[0])); sh = abs(float(qa[1] - qb[1]))
        checks.append({"metric": col, "small_ci": qa.tolist(), "large_ci": qb.tolist(), "large_width": width, "endpoint_tolerance": tol, "shift_low": sl, "shift_high": sh, "pass": bool(sl <= tol and sh <= tol)})
    passed = all(x["pass"] for x in checks)
    result = {"status": "PASS" if passed else "FAIL", "passed": passed, "n_small": n_small, "n_large": n_large, "checks": checks, "run_info": info}
    atomic_write_json(out / "manifests" / "BOOTSTRAP_CONVERGENCE_v321.json", result)
    st = _read_state(out); st["bootstrap_convergence"] = result; _write_state(out, st)
    if not passed:
        raise RuntimeError("v3.2.1 bootstrap convergence gate failed; do not run production")
    return result


def _require_convergence(out: Path) -> None:
    p = out / "manifests" / "BOOTSTRAP_CONVERGENCE_v321.json"
    if not p.exists() or not json.loads(p.read_text(encoding="utf-8")).get("passed", False):
        raise RuntimeError("v3.2.1 bootstrap convergence gate is not PASS")


def _report_markdown(summary: dict[str, pd.DataFrame], decisions: pd.DataFrame, calibration: dict[str, Any]) -> str:
    lines = [
        "# FGT correction-aware critical-scaling audit v3.2.1",
        "",
        "This report is post-processing only. It does not alter the accepted Monte Carlo trajectories.",
        "Tc used for nu inference is obtained from Binder and xi/L crossings and is not bounded by the legacy fixed-nu susceptibility-shift confidence interval.",
        "Quenched-realization percentile bootstrap intervals and finite-size/systematic robustness envelopes are reported separately.",
        "",
        f"Pristine calibration: **{calibration['status']}**",
        "",
        "## Decisions",
        "",
    ]
    for _, r in decisions.sort_values("p").iterrows():
        lines.append(f"- p={r.p:.2f}: **{r.decision}** — {r.reason}")
    lines += [
        "",
        "## Interpretation guard",
        "",
        "EVIDENCE_AGAINST_NU1, if it occurs, is evidence within the tested finite-size framework only. It is not proof of a new universality class. Marginal-disorder logarithmic corrections remain a competing explanation unless separately resolved with larger size range.",
    ]
    return "\n".join(lines) + "\n"


def run_real_audit(project_root: Path, spec_path: Path | None = None, workers: int = 4, force: bool = False) -> dict[str, Any]:
    c = ensure_central(project_root, spec_path)
    spec, out = c["spec"], c["out"]
    _require_synthetic(out); _require_convergence(out)
    nboot = int(spec["bootstrap"]["n_large"])
    namespace = f"production_v321_n{nboot}"
    info = run_bootstrap_chunks(project_root=project_root, out_dir=out, spec_path=spec_path, spec=spec, central=c["central"], n_bootstrap=nboot, workers=workers, cases=CASE_LABELS, namespace=namespace, force=force)
    boot = collect_bootstrap(out, checkpoint_root=Path(info["checkpoint_root"]), namespace=namespace, labels=CASE_LABELS, n_bootstrap=nboot)
    summary = summarize_production(boot, c["central"], spec)
    calibration = pristine_calibration(summary, spec, c["ctx"]["core_audit"])
    decisions = decide_cases(summary, spec, calibration)
    tables = out / "tables"
    atomic_write_csv(tables / "TABLE_TC_RG_CROSSING_BOOTSTRAP.csv", summary["tc"])
    atomic_write_csv(tables / "TABLE_NU_CORRECTION_AWARE.csv", summary["nu"])
    atomic_write_csv(tables / "TABLE_PAIRED_LMIN_DRIFT.csv", summary["drift"])
    atomic_write_csv(tables / "TABLE_EXPONENT_RATIO_CONSISTENCY.csv", summary["ratios"])
    atomic_write_csv(tables / "TABLE_CHANNEL_DELTA_NU.csv", summary["channel_delta"])
    atomic_write_csv(tables / "TABLE_FINAL_DECISIONS_v321.csv", decisions)
    atomic_write_csv(tables / "CLAIM_SCOPE_TABLE_v321.csv", claim_scope_table(decisions))
    atomic_write_json(out / "manifests" / "PRISTINE_CALIBRATION_v321.json", calibration)
    figcat = build_figures(summary, decisions, out / "figures", spec)
    atomic_write_csv(out / "figures" / "FIGURE_CATALOG_v321.csv", figcat)
    atomic_write_text(out / "FINAL_CORRECTION_AWARE_AUDIT_REPORT.md", _report_markdown(summary, decisions, calibration))
    result = {"status": "PASS", "created_utc": utc_now(), "n_bootstrap": nboot, "workers": int(workers), "pristine_calibration": calibration, "case_decisions": decisions.to_dict(orient="records"), "run_info": info}
    atomic_write_json(out / "FINAL_AUDIT_RESULT_v321.json", result)
    st = _read_state(out); st["production"] = result; _write_state(out, st)
    return result


def status(project_root: Path, spec_path: Path | None = None) -> dict[str, Any]:
    spec, out, _ = initialize_run(project_root, spec_path)
    cp = checkpoint_workspace(project_root, out, spec)
    progress: dict[str, Any] = {"output_dir": str(out), "checkpoint_root": str(cp), "state": _read_state(out), "namespaces": {}}
    if cp.exists():
        for ns in sorted(x for x in cp.iterdir() if x.is_dir()):
            progress["namespaces"][ns.name] = {d.name: len(list(d.glob("boot_*.json"))) for d in ns.iterdir() if d.is_dir()}
    return progress


def rebuild_figures(project_root: Path, spec_path: Path | None = None) -> dict[str, Any]:
    spec, out, _ = initialize_run(project_root, spec_path)
    tables = out / "tables"
    names = {
        "tc": "TABLE_TC_RG_CROSSING_BOOTSTRAP.csv", "nu": "TABLE_NU_CORRECTION_AWARE.csv",
        "drift": "TABLE_PAIRED_LMIN_DRIFT.csv", "ratios": "TABLE_EXPONENT_RATIO_CONSISTENCY.csv",
        "channel_delta": "TABLE_CHANNEL_DELTA_NU.csv",
    }
    summary = {k: pd.read_csv(tables / v) for k, v in names.items()}
    decisions = pd.read_csv(tables / "TABLE_FINAL_DECISIONS_v321.csv")
    cat = build_figures(summary, decisions, out / "figures", spec)
    atomic_write_csv(out / "figures" / "FIGURE_CATALOG_v321.csv", cat)
    return {"status": "PASS", "figure_files": int(len(cat))}


def validate_release(project_root: Path, spec_path: Path | None = None) -> dict[str, Any]:
    spec, out, _ = initialize_run(project_root, spec_path)
    required = [
        out / "manifests" / "MONTE_CARLO_CORE_AUDIT.json",
        out / "synthetic" / "SYNTHETIC_CHALLENGE_RESULT.json",
        out / "manifests" / "BOOTSTRAP_CONVERGENCE_v321.json",
        out / "manifests" / "PRISTINE_CALIBRATION_v321.json",
        out / "FINAL_AUDIT_RESULT_v321.json",
        out / "FINAL_CORRECTION_AWARE_AUDIT_REPORT.md",
    ]
    for name in ["TABLE_TC_RG_CROSSING_BOOTSTRAP.csv", "TABLE_NU_CORRECTION_AWARE.csv", "TABLE_PAIRED_LMIN_DRIFT.csv", "TABLE_EXPONENT_RATIO_CONSISTENCY.csv", "TABLE_CHANNEL_DELTA_NU.csv", "TABLE_FINAL_DECISIONS_v321.csv", "CLAIM_SCOPE_TABLE_v321.csv"]:
        required.append(out / "tables" / name)
    for stem in ["FIG01_rg_crossing_tc_vs_p", "FIG02_nu_bootstrap_and_robustness", "FIG03_paired_Lmin_drift", "FIG04_exponent_ratio_consistency", "FIG05_decision_summary"]:
        for fmt in spec["figures"]["formats"]:
            required.append(out / "figures" / f"{stem}.{fmt}")
    checks = [{"path": str(p), "pass": p.exists() and p.stat().st_size > 0} for p in required]
    # The production run is valid even if pristine calibration is FAIL; that is a scientific result.
    passed = all(x["pass"] for x in checks)
    atomic_write_csv(out / "RELEASE_VALIDATION_CHECKS_v321.csv", pd.DataFrame(checks))
    files = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and "checkpoints" not in p.parts and p.name != "OUTPUT_SHA256_MANIFEST_v321.csv":
            files.append({"path": p.relative_to(out).as_posix(), "sha256": sha256_file(p), "bytes": p.stat().st_size})
    atomic_write_csv(out / "OUTPUT_SHA256_MANIFEST_v321.csv", pd.DataFrame(files))
    result = {"status": "PASS" if passed else "FAIL", "passed": passed, "n_checks": len(checks), "n_failed": sum(not x["pass"] for x in checks)}
    atomic_write_json(out / "RELEASE_VALIDATION_RESULT_v321.json", result)
    if not passed:
        raise RuntimeError("v3.2.1 release validation failed")
    return result


def run_all(project_root: Path, spec_path: Path | None = None, workers: int = 4, force: bool = False) -> dict[str, Any]:
    syn = run_synthetic(project_root, spec_path, workers=workers)
    conv = run_bootstrap_convergence(project_root, spec_path, workers=workers, force=force)
    prod = run_real_audit(project_root, spec_path, workers=workers, force=force)
    rel = validate_release(project_root, spec_path)
    return {"synthetic": syn, "bootstrap_convergence": conv, "production": prod, "release_validation": rel}
