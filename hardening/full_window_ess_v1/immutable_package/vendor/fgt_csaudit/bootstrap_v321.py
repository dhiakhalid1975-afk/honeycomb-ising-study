from __future__ import annotations

from concurrent.futures import as_completed
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .analysis_v321 import PRIMARY
from .config import load_spec
from .correction_fit import exponent_ratio_at_tc, fit_nu_fixed_tc, restrict_curves
from .io import load_all_cases, validate_project_inputs
from .parallel import backend_name, make_executor
from .pb import SupportLock
from .provenance import atomic_write_json, safe_fsync, sha256_file, stable_json_hash
from .rg_tc import estimate_joint_tc
from .runtime import checkpoint_workspace

_WORKER_CASES = None
_WORKER_SPEC = None


def _worker_init(project_root: str, spec_path: str | None) -> None:
    global _WORKER_CASES, _WORKER_SPEC
    spec = load_spec(Path(spec_path) if spec_path else None)
    paths, _ = validate_project_inputs(Path(project_root), spec)
    _WORKER_CASES = load_all_cases(paths, spec)
    _WORKER_SPEC = spec


def _support_from_jsonable(d: dict[str, Any]) -> SupportLock:
    return SupportLock(
        sizes=tuple(int(x) for x in d["sizes"]),
        target_indices={int(k): np.asarray(v, dtype=int) for k, v in d["target_indices"].items()},
        x_low=float(d["x_low"]), x_high=float(d["x_high"]),
        n_ordered_residuals=int(d["n_ordered_residuals"]),
        reference_tc=float(d["reference_tc"]), reference_nu=float(d["reference_nu"]),
        x_window=tuple(float(x) for x in d["x_window"]), edge_points=int(d["edge_points"]),
    )


def _one_draw(case: Any, b: int, task: dict[str, Any], spec: dict[str, Any]) -> dict[str, Any]:
    curves = case.bootstrap_curves(int(b), int(spec["base_seed"]))
    branch_center = float(task["branch_center"])
    full_sizes = tuple(int(x) for x in spec["size_windows"]["full"])
    drop_sizes = tuple(int(x) for x in spec["size_windows"]["drop_smallest"])

    tc_full = estimate_joint_tc(curves, branch_center=branch_center, spec=spec, sizes=full_sizes)
    tc_drop = estimate_joint_tc(curves, branch_center=branch_center, spec=spec, sizes=drop_sizes)
    row: dict[str, Any] = {
        "case_label": case.label, "p": case.p, "bootstrap_index": int(b),
        "tc_full_joint": float(tc_full["joint_tc"]),
        "tc_full_binder": float(tc_full["binder_roa"]["estimate"]),
        "tc_full_xi": float(tc_full["xi_over_L"]["estimate"]),
        "tc_full_channel_spread": float(tc_full["channel_spread"]),
        "tc_drop_joint": float(tc_drop["joint_tc"]),
        "tc_drop_binder": float(tc_drop["binder_roa"]["estimate"]),
        "tc_drop_xi": float(tc_drop["xi_over_L"]["estimate"]),
        "tc_drop_channel_spread": float(tc_drop["channel_spread"]),
    }

    for wname, sizes, tcinfo in (("full", full_sizes, tc_full), ("drop", drop_sizes, tc_drop)):
        tc = float(tcinfo["joint_tc"])
        for ch in PRIMARY:
            key = f"{wname}:{ch}"
            support = _support_from_jsonable(task["supports"][key])
            bounds = tuple(float(x) for x in task["nu_bounds"][key])
            c = restrict_curves(curves[ch], sizes)
            if np.isfinite(tc):
                fr = fit_nu_fixed_tc(c, tc=tc, channel=ch, support=support, spec=spec, nu_bounds=bounds)
                row[f"nu_{ch}_{wname}"] = float(fr.nu)
                row[f"pb_{ch}_{wname}"] = float(fr.pb)
                row[f"boundary_{ch}_{wname}"] = bool(fr.boundary_hit)
                row[f"valid_{ch}_{wname}"] = bool(fr.valid)
            else:
                row[f"nu_{ch}_{wname}"] = float("nan")
                row[f"pb_{ch}_{wname}"] = float("nan")
                row[f"boundary_{ch}_{wname}"] = True
                row[f"valid_{ch}_{wname}"] = False

    # Tc-channel-choice sensitivity at full size window, using the same bootstrap draw.
    for ch in PRIMARY:
        support = _support_from_jsonable(task["supports"][f"full:{ch}"])
        bounds = tuple(float(x) for x in task["nu_bounds"][f"full:{ch}"])
        c = restrict_curves(curves[ch], full_sizes)
        for source, tc in (("binderTc", float(tc_full["binder_roa"]["estimate"])), ("xiTc", float(tc_full["xi_over_L"]["estimate"]))):
            if np.isfinite(tc):
                fr = fit_nu_fixed_tc(c, tc=tc, channel=ch, support=support, spec=spec, nu_bounds=bounds)
                row[f"nu_{ch}_{source}"] = float(fr.nu)
            else:
                row[f"nu_{ch}_{source}"] = float("nan")

    tc = float(tc_full["joint_tc"])
    for ch, kind in (("abs_m", "magnetization"), ("chi_abs", "susceptibility")):
        rr = exponent_ratio_at_tc(restrict_curves(curves[ch], full_sizes), tc, kind=kind) if np.isfinite(tc) else {"success": False, "ratio": float("nan")}
        row[f"ratio_{ch}_full"] = float(rr.get("ratio", float("nan")))
        rr2 = exponent_ratio_at_tc(restrict_curves(curves[ch], drop_sizes), float(tc_drop["joint_tc"]), kind=kind) if np.isfinite(float(tc_drop["joint_tc"])) else {"success": False, "ratio": float("nan")}
        row[f"ratio_{ch}_drop"] = float(rr2.get("ratio", float("nan")))
    return row


def _chunk_worker(task: dict[str, Any]) -> dict[str, Any]:
    if _WORKER_CASES is None or _WORKER_SPEC is None:
        raise RuntimeError("worker not initialized")
    case = _WORKER_CASES[task["case_label"]]
    rows = [_one_draw(case, b, task, _WORKER_SPEC) for b in range(int(task["start"]), int(task["stop"]))]
    return {"case_label": case.label, "p": case.p, "start": int(task["start"]), "stop": int(task["stop"]), "rows": rows}


def _chunk_paths(root: Path, label: str, start: int, stop: int) -> tuple[Path, Path]:
    d = root / label
    stem = f"boot_{start:05d}_{stop-1:05d}"
    return d / f"{stem}.csv", d / f"{stem}.json"


def _tmp(parent: Path, suffix: str) -> Path:
    parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".fg_", suffix=suffix, dir=str(parent))
    os.close(fd)
    return Path(name)


def _write_chunk(csv_path: Path, meta_path: Path, result: dict[str, Any], signature: str) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    t = _tmp(csv_path.parent, ".csv.tmp")
    try:
        pd.DataFrame(result["rows"]).to_csv(t, index=False)
        with t.open("rb") as fh:
            safe_fsync(fh)
        os.replace(t, csv_path)
    finally:
        t.unlink(missing_ok=True)
    atomic_write_json(meta_path, {
        "run_signature": signature, "case_label": result["case_label"],
        "start": result["start"], "stop": result["stop"],
        "csv_sha256": sha256_file(csv_path),
    })


def _valid(csv_path: Path, meta_path: Path, signature: str) -> bool:
    if not csv_path.exists() or not meta_path.exists():
        return False
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return meta.get("run_signature") == signature and meta.get("csv_sha256") == sha256_file(csv_path)
    except Exception:
        return False


def run_bootstrap_chunks(
    *, project_root: Path, out_dir: Path, spec_path: Path | None, spec: dict[str, Any],
    central: dict[str, Any], n_bootstrap: int, workers: int,
    cases: list[str], namespace: str, force: bool = False,
) -> dict[str, Any]:
    cp_root = checkpoint_workspace(project_root, out_dir, spec)
    chunk_root = cp_root / namespace
    chunk_root.mkdir(parents=True, exist_ok=True)
    print(f"[bootstrap-v3.2.1] checkpoint workspace: {chunk_root}", flush=True)
    _, hashes = validate_project_inputs(project_root, spec)
    signature = stable_json_hash({
        "spec": spec, "inputs": hashes, "n_bootstrap": int(n_bootstrap),
        "namespace": namespace, "backend": backend_name(),
        "method": "independent_rg_tc_then_scalar_nu_plus_paired_Lmin",
    })
    chunk_size = int(spec["bootstrap"]["chunk_size"])
    tasks: list[dict[str, Any]] = []
    skipped = 0
    for label in cases:
        cc = central["cases"][label]
        supports: dict[str, Any] = {}
        bounds: dict[str, Any] = {}
        for wold, wnew in (("full", "full"), ("drop_smallest", "drop")):
            for ch in PRIMARY:
                supports[f"{wnew}:{ch}"] = cc["supports"][wold][ch].to_jsonable()
                bounds[f"{wnew}:{ch}"] = list(cc["nu_bounds"][wold][ch])
        for start in range(0, int(n_bootstrap), chunk_size):
            stop = min(int(n_bootstrap), start + chunk_size)
            cp = _chunk_paths(chunk_root, label, start, stop)
            if not force and _valid(*cp, signature):
                skipped += 1
                continue
            tasks.append({
                "case_label": label, "start": start, "stop": stop,
                "branch_center": cc["branch_center"], "supports": supports, "nu_bounds": bounds,
            })
    failures: list[str] = []
    completed = 0
    if tasks:
        with make_executor(workers, initializer=_worker_init, initargs=(str(project_root), str(spec_path) if spec_path else None)) as ex:
            futs = {ex.submit(_chunk_worker, t): t for t in tasks}
            for fut in as_completed(futs):
                task = futs[fut]
                try:
                    res = fut.result()
                    cp = _chunk_paths(chunk_root, task["case_label"], task["start"], task["stop"])
                    _write_chunk(*cp, res, signature)
                    completed += 1
                    print(f"[bootstrap-v3.2.1] {task['case_label']} {task['start']}..{task['stop']-1} complete ({completed}/{len(tasks)})", flush=True)
                except Exception as exc:
                    failures.append(f"{task['case_label']}:{task['start']}-{task['stop']}: {type(exc).__name__}: {exc}")
                    print("[bootstrap-v3.2.1] FAILED " + failures[-1], flush=True)
    if failures:
        raise RuntimeError("bootstrap chunk failure:\n" + "\n".join(failures[:12]))
    return {
        "status": "PASS", "run_signature": signature, "checkpoint_root": str(cp_root),
        "chunk_root": str(chunk_root), "requested_tasks": len(tasks), "skipped_tasks": skipped,
        "completed_tasks": completed, "n_bootstrap": int(n_bootstrap), "workers": int(workers),
        "executor_backend": backend_name(),
    }


def collect_bootstrap(out_dir: Path, *, checkpoint_root: Path, namespace: str, labels: list[str], n_bootstrap: int) -> pd.DataFrame:
    root = checkpoint_root / namespace
    parts: list[pd.DataFrame] = []
    expected = set(range(int(n_bootstrap)))
    for label in labels:
        files = sorted((root / label).glob("boot_*.csv"))
        if not files:
            raise FileNotFoundError(f"no bootstrap chunks for {label}")
        df = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
        df = df.loc[df["bootstrap_index"] < int(n_bootstrap)].copy()
        got = set(int(x) for x in df["bootstrap_index"])
        if got != expected:
            raise RuntimeError(f"bootstrap coverage failure for {label}: missing={sorted(expected-got)[:10]}")
        parts.append(df)
    return pd.concat(parts, ignore_index=True)
