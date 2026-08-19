from __future__ import annotations

import argparse
from concurrent.futures import as_completed
import hashlib
import json
import math
import os
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

METHOD_VERSION = "NU_ASYMMETRIC_SENSITIVITY_v1_STRICT"
EXPECTED_SPEC_NAME = "FGT_CORRECTION_AWARE_CRITICAL_AUDIT_v3.2.1"
EXPECTED_SPEC_SHA256 = "98f5b9bee3b4f39495d70c88ed0b922f102e12abb12df0324c0975bffb3aace2"
PRIMARY = ("binder_roa", "xi_over_L")
CASE_LABELS = ("random_p080", "random_p085", "random_p090", "pristine_p100")
SCAN_POINTS = 181
BOUNDARY_FRACTION_OF_SEGMENT = 0.02
REPLAY_TOL = 5e-10

# Original-package functions/classes are installed at runtime from the exact package next to SPEC_LOCK.json.
_PB_RESIDUE = None
_BUILD_SUPPORT = None
_SUPPORT_LOCK = None
_ORIG_FIT = None
_RESTRICT = None
_EST_TC = None
_LOAD_SPEC = None
_VALIDATE_INPUTS = None
_LOAD_CASES = None
_CENTRAL_ALL = None
_MAKE_EXECUTOR = None
_BACKEND_NAME = None

# Worker state.
_W_CASES = None
_W_SPEC = None
_W_TASK_BUNDLE = None


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".nu_asym_", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    try:
        Path(tmp).write_text(text, encoding="utf-8")
        os.replace(tmp, path)
    finally:
        Path(tmp).unlink(missing_ok=True)


def atomic_json(path: Path, obj: Any) -> None:
    atomic_text(path, json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def atomic_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".nu_asym_", suffix=".csv.tmp", dir=str(path.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False)
        os.replace(tmp, path)
    finally:
        Path(tmp).unlink(missing_ok=True)


def is_within(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def install_original_modules(package_src: Path) -> None:
    global _PB_RESIDUE, _BUILD_SUPPORT, _SUPPORT_LOCK, _ORIG_FIT, _RESTRICT
    global _EST_TC, _LOAD_SPEC, _VALIDATE_INPUTS, _LOAD_CASES, _CENTRAL_ALL
    global _MAKE_EXECUTOR, _BACKEND_NAME

    src = package_src.resolve()
    if not (src / "fgt_csaudit").is_dir():
        raise RuntimeError(f"FAIL-CLOSED: package source has no fgt_csaudit directory: {src}")
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from fgt_csaudit.pb import SupportLock, build_locked_support, pb_residue
    from fgt_csaudit.correction_fit import fit_nu_fixed_tc, restrict_curves
    from fgt_csaudit.rg_tc import estimate_joint_tc
    from fgt_csaudit.config import load_spec
    from fgt_csaudit.io import load_all_cases, validate_project_inputs
    from fgt_csaudit.analysis_v321 import central_all
    from fgt_csaudit.parallel import backend_name, make_executor

    _PB_RESIDUE = pb_residue
    _BUILD_SUPPORT = build_locked_support
    _SUPPORT_LOCK = SupportLock
    _ORIG_FIT = fit_nu_fixed_tc
    _RESTRICT = restrict_curves
    _EST_TC = estimate_joint_tc
    _LOAD_SPEC = load_spec
    _VALIDATE_INPUTS = validate_project_inputs
    _LOAD_CASES = load_all_cases
    _CENTRAL_ALL = central_all
    _MAKE_EXECUTOR = make_executor
    _BACKEND_NAME = backend_name


def support_from_jsonable(d: dict[str, Any]):
    return _SUPPORT_LOCK(
        sizes=tuple(int(x) for x in d["sizes"]),
        target_indices={int(k): np.asarray(v, dtype=int) for k, v in d["target_indices"].items()},
        x_low=float(d["x_low"]),
        x_high=float(d["x_high"]),
        n_ordered_residuals=int(d["n_ordered_residuals"]),
        reference_tc=float(d["reference_tc"]),
        reference_nu=float(d["reference_nu"]),
        x_window=tuple(float(x) for x in d["x_window"]),
        edge_points=int(d["edge_points"]),
    )


def segments_from_mask(grid: np.ndarray, mask: np.ndarray) -> list[tuple[float, float]]:
    grid = np.asarray(grid, float)
    mask = np.asarray(mask, bool)
    if len(grid) != len(mask):
        raise ValueError("grid/mask length mismatch")
    out: list[tuple[float, float]] = []
    start = None
    for i, ok in enumerate(mask):
        if ok and start is None:
            start = i
        end_run = start is not None and (not ok or i == len(mask) - 1)
        if end_run:
            end = i if ok and i == len(mask) - 1 else i - 1
            out.append((float(grid[start]), float(grid[end])))
            start = None
    return out


def scan_feasible_segments(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    channel: str,
    support: Any,
    spec: dict[str, Any],
    interpolation: str | None = None,
) -> dict[str, Any]:
    declared = tuple(float(x) for x in spec["nu_fit"]["declared_bounds"])
    method = str(interpolation or spec["nu_fit"]["primary_interpolation"])
    penalty = float(spec["nu_fit"]["invalid_penalty"])
    grid = np.linspace(declared[0], declared[1], SCAN_POINTS)
    vals = np.full(len(grid), np.nan, dtype=float)
    valid = np.zeros(len(grid), dtype=bool)
    for i, nu in enumerate(grid):
        val, _, ok = _PB_RESIDUE(
            curves,
            tc=float(tc),
            nu=float(nu),
            q=0.0,
            channel=channel,
            support=support,
            interpolation=method,
            invalid_penalty=penalty,
        )
        vals[i] = float(val)
        valid[i] = bool(ok and np.isfinite(val) and float(val) < 0.1 * penalty)
    segs = segments_from_mask(grid, valid)
    if not segs:
        raise RuntimeError(f"{channel}: no feasible nu segment in declared range {declared}")
    return {
        "declared_bounds": [declared[0], declared[1]],
        "scan_points": SCAN_POINTS,
        "segments": [[a, b] for a, b in segs],
        "n_segments": len(segs),
        "feasible_grid_count": int(np.sum(valid)),
        "feasible_min": float(min(a for a, _ in segs)),
        "feasible_max": float(max(b for _, b in segs)),
        "includes_nu1": bool(any(a <= 1.0 <= b for a, b in segs)),
    }


@dataclass(frozen=True)
class AsymFit:
    channel: str
    tc: float
    nu: float
    pb: float
    valid: bool
    boundary_hit: bool
    segment_index: int
    segment_bounds: tuple[float, float]
    n_segments: int
    nfev: int
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _segment_index(nu: float, segments: list[tuple[float, float]], tol: float = 1e-12) -> int:
    for i, (a, b) in enumerate(segments):
        if a - tol <= nu <= b + tol:
            return i
    return -1


def fit_nu_asymmetric(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    channel: str,
    support: Any,
    spec: dict[str, Any],
    segments: list[tuple[float, float]] | list[list[float]],
    interpolation: str | None = None,
) -> AsymFit:
    segs = [(float(a), float(b)) for a, b in segments]
    if not segs:
        return AsymFit(channel, float(tc), float("nan"), float("nan"), False, True, -1, (float("nan"), float("nan")), 0, 0, "no locked feasible segments")
    for a, b in segs:
        if a > b:
            raise ValueError("segment lower bound exceeds upper bound")
    method = str(interpolation or spec["nu_fit"]["primary_interpolation"])
    penalty = float(spec["nu_fit"]["invalid_penalty"])
    ncoarse = int(spec["nu_fit"]["coarse_points"])
    hull_lo = min(a for a, _ in segs)
    hull_hi = max(b for _, b in segs)
    if not hull_lo < hull_hi:
        return AsymFit(channel, float(tc), float(hull_lo), float(penalty), False, True, 0, segs[0], len(segs), 0, "degenerate feasible hull")

    def f(nu: float) -> float:
        return float(_PB_RESIDUE(
            curves,
            tc=float(tc),
            nu=float(nu),
            q=0.0,
            channel=channel,
            support=support,
            interpolation=method,
            invalid_penalty=penalty,
        )[0])

    grid = np.linspace(hull_lo, hull_hi, ncoarse)
    allowed = np.asarray([_segment_index(float(x), segs) >= 0 for x in grid], dtype=bool)
    vals = np.full(len(grid), float(penalty), dtype=float)
    nfev = 0
    for i in np.flatnonzero(allowed):
        vals[i] = f(float(grid[i]))
        nfev += 1
    valid = allowed & np.isfinite(vals) & (vals < 0.1 * penalty)
    if not np.any(valid):
        return AsymFit(channel, float(tc), float("nan"), float(penalty), False, True, -1, (hull_lo, hull_hi), len(segs), nfev, "no valid coarse points inside locked asymmetric segments")

    idx_valid = np.flatnonzero(valid)
    ibest = int(idx_valid[np.argmin(vals[idx_valid])])
    best_nu = float(grid[ibest])
    best_pb = float(vals[ibest])
    sidx = _segment_index(best_nu, segs)
    if sidx < 0:
        raise RuntimeError("software integrity failure: best coarse point is outside feasible segments")
    slo, shi = segs[sidx]

    left_i = max(0, ibest - 1)
    right_i = min(len(grid) - 1, ibest + 1)
    a = max(slo, float(grid[left_i]))
    b = min(shi, float(grid[right_i]))
    if a < b:
        try:
            res = minimize_scalar(
                f,
                bounds=(a, b),
                method="bounded",
                options={"xatol": float(spec["nu_fit"]["refine_xatol"])},
            )
            nfev += int(getattr(res, "nfev", 0))
            score = float(f(float(res.x)))
            nfev += 1
            if bool(res.success) and np.isfinite(score) and score < 0.1 * penalty and score < best_pb:
                best_nu, best_pb = float(res.x), score
        except Exception:
            pass

    width = shi - slo
    if width <= 0:
        boundary = True
    else:
        seg_allowed_idx = np.asarray([i for i, x in enumerate(grid) if slo - 1e-12 <= x <= shi + 1e-12], dtype=int)
        coarse_edge = bool(len(seg_allowed_idx) and ibest in {int(seg_allowed_idx[0]), int(seg_allowed_idx[-1])})
        boundary = bool(
            (best_nu - slo) <= BOUNDARY_FRACTION_OF_SEGMENT * width
            or (shi - best_nu) <= BOUNDARY_FRACTION_OF_SEGMENT * width
            or coarse_edge
        )
    return AsymFit(
        channel=channel,
        tc=float(tc),
        nu=float(best_nu),
        pb=float(best_pb),
        valid=True,
        boundary_hit=boundary,
        segment_index=int(sidx),
        segment_bounds=(float(slo), float(shi)),
        n_segments=len(segs),
        nfev=int(nfev),
        message="locked asymmetric feasible segments + deterministic coarse grid + bounded scalar refinement",
    )


def synthetic_curves(tc: float, nu: float) -> dict[str, dict[int, tuple[np.ndarray, np.ndarray]]]:
    sizes = (40, 60, 80, 100, 120)
    out = {"binder_roa": {}, "xi_over_L": {}}
    for L in sizes:
        T = np.linspace(tc - 0.075, tc + 0.075, 31)
        x = (T - tc) / tc * L ** (1.0 / nu)
        out["binder_roa"][L] = (T, 0.610 - 0.105 * np.tanh(x / 2.2))
        out["xi_over_L"][L] = (T, 0.920 - 0.170 * np.tanh(x / 2.0))
    return out


def run_selftests(spec: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    grid = np.linspace(0.0, 0.4, 5)
    seg = segments_from_mask(grid, np.asarray([True, True, False, True, True]))
    checks.append({"check": "topology_two_segments", "pass": seg == [(0.0, 0.1), (0.30000000000000004, 0.4)], "detail": seg})

    for true_nu in (1.0, 1.12):
        curves = synthetic_curves(tc=1.2, nu=true_nu)
        for ch in PRIMARY:
            support = _BUILD_SUPPORT(
                curves[ch],
                tc=1.2,
                nu=1.0,
                channel=ch,
                q=0.0,
                x_window=tuple(float(x) for x in spec["support"]["x_window"]),
                edge_points=int(spec["support"]["edge_points"]),
                minimum_points_per_size=int(spec["support"]["minimum_points_per_size"]),
            )
            scan = scan_feasible_segments(curves[ch], tc=1.2, channel=ch, support=support, spec=spec)
            fit = fit_nu_asymmetric(curves[ch], tc=1.2, channel=ch, support=support, spec=spec, segments=scan["segments"])
            tol = 0.035 if true_nu == 1.0 else 0.05
            ok = bool(fit.valid and np.isfinite(fit.nu) and abs(fit.nu - true_nu) <= tol)
            checks.append({
                "check": f"recover_nu_{true_nu:.2f}_{ch}",
                "pass": ok,
                "detail": {"true_nu": true_nu, "fit_nu": fit.nu, "segments": scan["segments"], "boundary": fit.boundary_hit},
            })

    passed = all(bool(x["pass"]) for x in checks)
    return {"status": "PASS" if passed else "FAIL", "passed": passed, "n_checks": len(checks), "checks": checks}


def snapshot(paths: list[Path]) -> dict[str, Any]:
    rec: dict[str, Any] = {}
    for p in sorted({Path(x).resolve() for x in paths}, key=lambda x: str(x).lower()):
        if not p.exists() or not p.is_file():
            raise FileNotFoundError(f"snapshot path missing: {p}")
        st = p.stat()
        rec[str(p)] = {"sha256": sha256_file(p), "bytes": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}
    return rec


def compare_snapshots(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    keys = sorted(set(before) | set(after))
    rows = []
    for k in keys:
        b = before.get(k)
        a = after.get(k)
        same = bool(b == a)
        rows.append({"path": k, "pass": same, "before": b, "after": a})
    passed = all(r["pass"] for r in rows)
    return {"status": "PASS" if passed else "FAIL", "passed": passed, "checks": rows}


def required_result_paths(results: Path) -> list[Path]:
    req = [
        results / "tables" / "TABLE_TC_RG_CROSSING_BOOTSTRAP.csv",
        results / "tables" / "TABLE_NU_CORRECTION_AWARE.csv",
        results / "tables" / "TABLE_PAIRED_LMIN_DRIFT.csv",
        results / "tables" / "TABLE_FINAL_DECISIONS_v321.csv",
        results / "tables" / "LOCKED_SUPPORT_AND_SYMMETRIC_NU_BOUNDS.json",
        results / "FINAL_AUDIT_RESULT_v321.json",
        results / "manifests" / "PRISTINE_CALIBRATION_v321.json",
    ]
    missing = [str(p) for p in req if not p.exists()]
    if missing:
        raise FileNotFoundError("required original result files are missing:\n" + "\n".join(missing))
    return req


def critical_source_paths(package_root: Path) -> list[Path]:
    base = package_root / "src" / "fgt_csaudit"
    names = [
        "pb.py", "correction_fit.py", "rg_tc.py", "analysis_v321.py", "bootstrap_v321.py",
        "io.py", "config.py", "parallel.py", "decision_v321.py", "summarize_v321.py",
    ]
    req = [base / n for n in names]
    missing = [str(p) for p in req if not p.exists()]
    if missing:
        raise FileNotFoundError("required original source files are missing:\n" + "\n".join(missing))
    return req


def build_domain_map(central: dict[str, Any], spec: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    bundle: dict[str, Any] = {}
    windows = {
        "full": tuple(int(x) for x in spec["size_windows"]["full"]),
        "drop_smallest": tuple(int(x) for x in spec["size_windows"]["drop_smallest"]),
        "diagnostic_large_only": tuple(int(x) for x in spec["size_windows"]["diagnostic_large_only"]),
    }
    for label in CASE_LABELS:
        cc = central["cases"][label]
        bundle[label] = {
            "p": float(cc["p"]),
            "branch_center": float(cc["branch_center"]),
            "supports": {},
            "sym_bounds": {},
            "asym_segments": {},
        }
        for wname, sizes in windows.items():
            tcinfo = cc["tc_by_window"].get(wname)
            if not tcinfo or not np.isfinite(float(tcinfo["joint_tc"])):
                raise RuntimeError(f"{label}/{wname}: central Tc unavailable")
            tc = float(tcinfo["joint_tc"])
            for ch in PRIMARY:
                curves = _RESTRICT(cc["central_curves"][ch], sizes)
                support = cc["supports"][wname][ch]
                scan = scan_feasible_segments(curves, tc=tc, channel=ch, support=support, spec=spec)
                sym = tuple(float(x) for x in cc["nu_bounds"][wname][ch])
                key = f"{wname}:{ch}"
                bundle[label]["supports"][key] = support.to_jsonable()
                bundle[label]["sym_bounds"][key] = [sym[0], sym[1]]
                bundle[label]["asym_segments"][key] = scan["segments"]
                rows.append({
                    "case_label": label,
                    "p": float(cc["p"]),
                    "window": wname,
                    "channel": ch,
                    "tc_joint": tc,
                    "declared_lo": float(scan["declared_bounds"][0]),
                    "declared_hi": float(scan["declared_bounds"][1]),
                    "symmetric_lo": sym[0],
                    "symmetric_hi": sym[1],
                    "asymmetric_feasible_min": float(scan["feasible_min"]),
                    "asymmetric_feasible_max": float(scan["feasible_max"]),
                    "n_feasible_segments": int(scan["n_segments"]),
                    "segments_json": json.dumps(scan["segments"], separators=(",", ":")),
                    "feasible_grid_count": int(scan["feasible_grid_count"]),
                    "includes_nu1": bool(scan["includes_nu1"]),
                    "scan_points": SCAN_POINTS,
                })
    return pd.DataFrame(rows), bundle


def worker_init(project_root: str, spec_path: str, package_src: str, task_bundle: dict[str, Any]) -> None:
    global _W_CASES, _W_SPEC, _W_TASK_BUNDLE
    install_original_modules(Path(package_src))
    _W_SPEC = _LOAD_SPEC(Path(spec_path))
    paths, _ = _VALIDATE_INPUTS(Path(project_root), _W_SPEC)
    _W_CASES = _LOAD_CASES(paths, _W_SPEC)
    _W_TASK_BUNDLE = task_bundle


def one_draw(label: str, b: int) -> dict[str, Any]:
    case = _W_CASES[label]
    spec = _W_SPEC
    tb = _W_TASK_BUNDLE[label]
    curves = case.bootstrap_curves(int(b), int(spec["base_seed"]))
    full_sizes = tuple(int(x) for x in spec["size_windows"]["full"])
    drop_sizes = tuple(int(x) for x in spec["size_windows"]["drop_smallest"])
    tc_full = _EST_TC(curves, branch_center=float(tb["branch_center"]), spec=spec, sizes=full_sizes)
    tc_drop = _EST_TC(curves, branch_center=float(tb["branch_center"]), spec=spec, sizes=drop_sizes)
    row: dict[str, Any] = {
        "case_label": label,
        "p": float(case.p),
        "bootstrap_index": int(b),
        "tc_full_joint": float(tc_full["joint_tc"]),
        "tc_full_binder": float(tc_full["binder_roa"]["estimate"]),
        "tc_full_xi": float(tc_full["xi_over_L"]["estimate"]),
        "tc_drop_joint": float(tc_drop["joint_tc"]),
    }

    for wname, sizes, tcinfo, wkey in (
        ("full", full_sizes, tc_full, "full"),
        ("drop", drop_sizes, tc_drop, "drop_smallest"),
    ):
        tc = float(tcinfo["joint_tc"])
        for ch in PRIMARY:
            key = f"{wkey}:{ch}"
            support = support_from_jsonable(tb["supports"][key])
            c = _RESTRICT(curves[ch], sizes)
            sym_bounds = tuple(float(x) for x in tb["sym_bounds"][key])
            asym_segments = [(float(a), float(z)) for a, z in tb["asym_segments"][key]]
            if np.isfinite(tc):
                sf = _ORIG_FIT(c, tc=tc, channel=ch, support=support, spec=spec, nu_bounds=sym_bounds)
                af = fit_nu_asymmetric(c, tc=tc, channel=ch, support=support, spec=spec, segments=asym_segments)
                row[f"nu_sym_{ch}_{wname}"] = float(sf.nu)
                row[f"pb_sym_{ch}_{wname}"] = float(sf.pb)
                row[f"boundary_sym_{ch}_{wname}"] = bool(sf.boundary_hit)
                row[f"valid_sym_{ch}_{wname}"] = bool(sf.valid)
                row[f"nu_asym_{ch}_{wname}"] = float(af.nu)
                row[f"pb_asym_{ch}_{wname}"] = float(af.pb)
                row[f"boundary_asym_{ch}_{wname}"] = bool(af.boundary_hit)
                row[f"valid_asym_{ch}_{wname}"] = bool(af.valid)
                row[f"segment_asym_{ch}_{wname}"] = int(af.segment_index)
            else:
                for mode in ("sym", "asym"):
                    row[f"nu_{mode}_{ch}_{wname}"] = float("nan")
                    row[f"pb_{mode}_{ch}_{wname}"] = float("nan")
                    row[f"boundary_{mode}_{ch}_{wname}"] = True
                    row[f"valid_{mode}_{ch}_{wname}"] = False
                row[f"segment_asym_{ch}_{wname}"] = -1
    return row


def chunk_worker(task: dict[str, Any]) -> dict[str, Any]:
    label = str(task["case_label"])
    start = int(task["start"])
    stop = int(task["stop"])
    return {"case_label": label, "start": start, "stop": stop, "rows": [one_draw(label, b) for b in range(start, stop)]}


def chunk_paths(root: Path, label: str, start: int, stop: int) -> tuple[Path, Path]:
    d = root / label
    stem = f"boot_{start:05d}_{stop - 1:05d}"
    return d / f"{stem}.csv", d / f"{stem}.json"


def valid_checkpoint(csv_path: Path, meta_path: Path, signature: str) -> bool:
    if not csv_path.exists() or not meta_path.exists():
        return False
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return meta.get("run_signature") == signature and meta.get("csv_sha256") == sha256_file(csv_path)
    except Exception:
        return False


def write_checkpoint(csv_path: Path, meta_path: Path, result: dict[str, Any], signature: str) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_csv(csv_path, pd.DataFrame(result["rows"]))
    atomic_json(meta_path, {
        "run_signature": signature,
        "case_label": result["case_label"],
        "start": int(result["start"]),
        "stop": int(result["stop"]),
        "csv_sha256": sha256_file(csv_path),
    })


def run_bootstrap(
    *,
    project_root: Path,
    spec_path: Path,
    package_src: Path,
    spec: dict[str, Any],
    task_bundle: dict[str, Any],
    out: Path,
    workers: int,
    force: bool,
    signature: str,
) -> pd.DataFrame:
    nboot = int(spec["bootstrap"]["n_large"])
    chunk_size = int(spec["bootstrap"]["chunk_size"])
    cp_root = out / "checkpoints"
    tasks = []
    skipped = 0
    for label in CASE_LABELS:
        for start in range(0, nboot, chunk_size):
            stop = min(nboot, start + chunk_size)
            csvp, metap = chunk_paths(cp_root, label, start, stop)
            if not force and valid_checkpoint(csvp, metap, signature):
                skipped += 1
                continue
            tasks.append({"case_label": label, "start": start, "stop": stop})

    print(f"[bootstrap] backend={_BACKEND_NAME()} workers={workers} nboot={nboot} tasks={len(tasks)} skipped={skipped}", flush=True)
    failures: list[str] = []
    if tasks:
        with _MAKE_EXECUTOR(
            workers,
            initializer=worker_init,
            initargs=(str(project_root), str(spec_path), str(package_src), task_bundle),
        ) as ex:
            futs = {ex.submit(chunk_worker, t): t for t in tasks}
            done = 0
            for fut in as_completed(futs):
                task = futs[fut]
                try:
                    result = fut.result()
                    csvp, metap = chunk_paths(cp_root, task["case_label"], task["start"], task["stop"])
                    write_checkpoint(csvp, metap, result, signature)
                    done += 1
                    print(f"[bootstrap] {task['case_label']} {task['start']}..{task['stop']-1} complete ({done}/{len(tasks)})", flush=True)
                except Exception as exc:
                    failures.append(f"{task['case_label']}:{task['start']}-{task['stop']}: {type(exc).__name__}: {exc}")
    if failures:
        raise RuntimeError("bootstrap chunk failure:\n" + "\n".join(failures[:20]))

    parts: list[pd.DataFrame] = []
    expected = set(range(nboot))
    for label in CASE_LABELS:
        files = sorted((cp_root / label).glob("boot_*.csv"))
        if not files:
            raise FileNotFoundError(f"no bootstrap checkpoints for {label}")
        df = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
        df = df.loc[df["bootstrap_index"] < nboot].copy()
        got = set(int(x) for x in df["bootstrap_index"])
        if got != expected:
            raise RuntimeError(f"bootstrap coverage failure for {label}: missing={sorted(expected - got)[:10]}")
        if len(df) != nboot:
            raise RuntimeError(f"bootstrap duplicate-row failure for {label}: got {len(df)} rows, expected {nboot}")
        parts.append(df)
    all_df = pd.concat(parts, ignore_index=True).sort_values(["case_label", "bootstrap_index"]).reset_index(drop=True)
    atomic_csv(out / "tables" / "BOOTSTRAP_DRAW_LEVEL_SYM_VS_ASYM.csv", all_df)
    return all_df


def numeric_summary(values: Any) -> dict[str, float | int]:
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return {"n": 0, "median": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    return {
        "n": int(len(a)),
        "median": float(np.median(a)),
        "ci_low": float(np.quantile(a, 0.025)),
        "ci_high": float(np.quantile(a, 0.975)),
    }


def bool_array(s: pd.Series) -> np.ndarray:
    if s.dtype == bool:
        return s.to_numpy(bool)
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"]).to_numpy(bool)


def summarize_draws(boot: pd.DataFrame, domain_map: pd.DataFrame, spec: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    max_boundary = float(spec["decision"].get("max_primary_boundary_fraction", 0.10))
    min_valid = float(spec["decision"].get("min_primary_valid_fraction", 0.95))
    for label in CASE_LABELS:
        g = boot.loc[boot.case_label == label].copy()
        p = float(g.p.iloc[0])
        for ch in PRIMARY:
            ss = numeric_summary(g[f"nu_sym_{ch}_full"])
            aa = numeric_summary(g[f"nu_asym_{ch}_full"])
            ad = numeric_summary(g[f"nu_asym_{ch}_drop"])
            drift = g[f"nu_asym_{ch}_drop"].to_numpy(float) - g[f"nu_asym_{ch}_full"].to_numpy(float)
            ds = numeric_summary(drift)
            delta = g[f"nu_asym_{ch}_full"].to_numpy(float) - g[f"nu_sym_{ch}_full"].to_numpy(float)
            dx = numeric_summary(delta)
            sym_b = float(np.mean(bool_array(g[f"boundary_sym_{ch}_full"])))
            sym_v = float(np.mean(bool_array(g[f"valid_sym_{ch}_full"])))
            asym_b = float(np.mean(bool_array(g[f"boundary_asym_{ch}_full"])))
            asym_v = float(np.mean(bool_array(g[f"valid_asym_{ch}_full"])))
            toward = False
            significant = False
            med = float(aa["median"])
            if np.isfinite(med) and np.isfinite(float(ds["ci_low"])) and np.isfinite(float(ds["ci_high"])):
                if med > 1.0:
                    toward = float(ds["median"]) < 0.0
                    significant = float(ds["ci_high"]) < 0.0
                elif med < 1.0:
                    toward = float(ds["median"]) > 0.0
                    significant = float(ds["ci_low"]) > 0.0
            dm = domain_map.loc[(domain_map.case_label == label) & (domain_map.window == "full") & (domain_map.channel == ch)]
            if len(dm) != 1:
                raise RuntimeError(f"expected one domain row for {label}/{ch}, got {len(dm)}")
            d = dm.iloc[0]
            rows.append({
                "case_label": label,
                "p": p,
                "channel": ch,
                "n_bootstrap": int(aa["n"]),
                "symmetric_nu_median": float(ss["median"]),
                "symmetric_nu_ci_low": float(ss["ci_low"]),
                "symmetric_nu_ci_high": float(ss["ci_high"]),
                "asymmetric_nu_median": float(aa["median"]),
                "asymmetric_nu_ci_low": float(aa["ci_low"]),
                "asymmetric_nu_ci_high": float(aa["ci_high"]),
                "asymmetric_drop_nu_median": float(ad["median"]),
                "asymmetric_drop_nu_ci_low": float(ad["ci_low"]),
                "asymmetric_drop_nu_ci_high": float(ad["ci_high"]),
                "paired_asym_minus_sym_median": float(dx["median"]),
                "paired_asym_minus_sym_ci_low": float(dx["ci_low"]),
                "paired_asym_minus_sym_ci_high": float(dx["ci_high"]),
                "asymmetric_Lmin_drift_median": float(ds["median"]),
                "asymmetric_Lmin_drift_ci_low": float(ds["ci_low"]),
                "asymmetric_Lmin_drift_ci_high": float(ds["ci_high"]),
                "asymmetric_Lmin_drift_toward_one": bool(toward),
                "asymmetric_Lmin_drift_significant_toward_one": bool(significant),
                "symmetric_boundary_fraction": sym_b,
                "symmetric_valid_fraction": sym_v,
                "asymmetric_boundary_fraction": asym_b,
                "asymmetric_valid_fraction": asym_v,
                "asymmetric_identifiability_gate_pass": bool(asym_b <= max_boundary and asym_v >= min_valid),
                "asymmetric_ci_contains_1": bool(np.isfinite(float(aa["ci_low"])) and float(aa["ci_low"]) <= 1.0 <= float(aa["ci_high"])),
                "symmetric_feasible_lo": float(d.symmetric_lo),
                "symmetric_feasible_hi": float(d.symmetric_hi),
                "asymmetric_feasible_min": float(d.asymmetric_feasible_min),
                "asymmetric_feasible_max": float(d.asymmetric_feasible_max),
                "n_feasible_segments": int(d.n_feasible_segments),
                "segments_json": str(d.segments_json),
                "max_boundary_gate": max_boundary,
                "min_valid_gate": min_valid,
            })
    return pd.DataFrame(rows)


def replay_audit(boot: pd.DataFrame, summary: pd.DataFrame, original_results: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    orig_nu = pd.read_csv(original_results / "tables" / "TABLE_NU_CORRECTION_AWARE.csv")
    orig_tc = pd.read_csv(original_results / "tables" / "TABLE_TC_RG_CROSSING_BOOTSTRAP.csv")
    rows: list[dict[str, Any]] = []

    def add(kind: str, label: str, ch: str, metric: str, replay: float, original: float) -> None:
        diff = abs(float(replay) - float(original)) if np.isfinite(replay) and np.isfinite(original) else float("inf")
        rows.append({
            "kind": kind,
            "case_label": label,
            "channel": ch,
            "metric": metric,
            "replay": replay,
            "original": original,
            "abs_diff": diff,
            "tolerance": REPLAY_TOL,
            "pass": bool(diff <= REPLAY_TOL),
        })

    for label in CASE_LABELS:
        g = boot.loc[boot.case_label == label]
        ot = orig_tc.loc[orig_tc.case_label == label]
        if len(ot) != 1:
            raise RuntimeError(f"original Tc table row count is not one for {label}")
        ot = ot.iloc[0]
        ts = numeric_summary(g.tc_full_joint)
        add("tc", label, "joint", "tc_bootstrap_median", float(ts["median"]), float(ot.tc_bootstrap_median))
        add("tc", label, "joint", "tc_bootstrap_ci_low", float(ts["ci_low"]), float(ot.tc_bootstrap_ci_low))
        add("tc", label, "joint", "tc_bootstrap_ci_high", float(ts["ci_high"]), float(ot.tc_bootstrap_ci_high))

        for ch in PRIMARY:
            s = summary.loc[(summary.case_label == label) & (summary.channel == ch)]
            o = orig_nu.loc[(orig_nu.case_label == label) & (orig_nu.channel == ch)]
            if len(s) != 1 or len(o) != 1:
                raise RuntimeError(f"replay/original nu row count failure for {label}/{ch}: {len(s)}/{len(o)}")
            s = s.iloc[0]
            o = o.iloc[0]
            for sm, om in (
                ("symmetric_nu_median", "nu_bootstrap_median"),
                ("symmetric_nu_ci_low", "nu_bootstrap_ci_low"),
                ("symmetric_nu_ci_high", "nu_bootstrap_ci_high"),
                ("symmetric_boundary_fraction", "bootstrap_boundary_fraction"),
                ("symmetric_valid_fraction", "bootstrap_valid_fraction"),
            ):
                add("nu", label, ch, om, float(s[sm]), float(o[om]))

    df = pd.DataFrame(rows)
    passed = bool(df["pass"].all())
    result = {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "tolerance": REPLAY_TOL,
        "n_checks": int(len(df)),
        "n_failed": int((~df["pass"]).sum()),
        "meaning": "The symmetric replay must reproduce the original v3.2.1 Tc/nu bootstrap summaries before asymmetric sensitivity is interpreted.",
    }
    return df, result


def shadow_decisions(summary: pd.DataFrame, spec: dict[str, Any], replay_pass: bool) -> pd.DataFrame:
    rows = []
    pristine = str(spec["calibration"]["pristine_case"])
    for label in CASE_LABELS:
        g = summary.loc[summary.case_label == label]
        if len(g) != 2:
            raise RuntimeError(f"expected two primary-channel summary rows for {label}")
        gates = [bool(x) for x in g.asymmetric_identifiability_gate_pass]
        drift = [bool(x) for x in g.asymmetric_Lmin_drift_significant_toward_one]
        ci1 = [bool(x) for x in g.asymmetric_ci_contains_1]
        if not replay_pass:
            status = "INDETERMINATE_REPLAY_MISMATCH"
            reason = "symmetric replay did not reproduce original locked production summaries"
        elif label == pristine:
            if all(gates) and all(ci1):
                status = "PRISTINE_DIAGNOSTIC_PASS"
                reason = "both asymmetric primary gates pass and both nominal bootstrap intervals contain nu=1; calibration is not redefined here"
            else:
                status = "PRISTINE_DIAGNOSTIC_REVIEW"
                reason = "asymmetric sensitivity changes a pristine diagnostic; do not alter original calibration automatically"
        elif any(drift):
            status = "ROBUST_NONIDENTIFIABLE_CORRECTION_DOMINATED"
            reason = "at least one primary channel retains significant paired Lmin drift toward nu=1 under asymmetric domains"
        elif not all(gates):
            status = "ROBUST_NONIDENTIFIABLE_LIMITED_RANGE"
            reason = "at least one primary channel still fails the locked boundary/validity identifiability gate under asymmetric domains"
        else:
            status = "REQUIRES_FULL_DECISION_REVIEW"
            reason = "both asymmetric primary gates pass and no significant paired drift veto remains; this sensitivity alone cannot replace the original decision layer"
        rows.append({
            "case_label": label,
            "p": float(g.p.iloc[0]),
            "shadow_status": status,
            "reason": reason,
            "binder_gate_pass": bool(g.loc[g.channel == "binder_roa", "asymmetric_identifiability_gate_pass"].iloc[0]),
            "xi_gate_pass": bool(g.loc[g.channel == "xi_over_L", "asymmetric_identifiability_gate_pass"].iloc[0]),
            "binder_boundary_fraction": float(g.loc[g.channel == "binder_roa", "asymmetric_boundary_fraction"].iloc[0]),
            "xi_boundary_fraction": float(g.loc[g.channel == "xi_over_L", "asymmetric_boundary_fraction"].iloc[0]),
            "scope": "shadow sensitivity only; never overwrites TABLE_FINAL_DECISIONS_v321.csv",
        })
    return pd.DataFrame(rows)


def make_figure(summary: pd.DataFrame, out: Path, spec: dict[str, Any]) -> dict[str, Any]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    max_boundary = float(spec["decision"].get("max_primary_boundary_fraction", 0.10))
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.2), sharex="col")
    labels = ["(a)", "(b)", "(c)", "(d)"]
    for j, ch in enumerate(PRIMARY):
        g = summary.loc[summary.channel == ch].sort_values("p")
        p = g.p.to_numpy(float)
        sm = g.symmetric_nu_median.to_numpy(float)
        sl = g.symmetric_nu_ci_low.to_numpy(float)
        sh = g.symmetric_nu_ci_high.to_numpy(float)
        am = g.asymmetric_nu_median.to_numpy(float)
        al = g.asymmetric_nu_ci_low.to_numpy(float)
        ah = g.asymmetric_nu_ci_high.to_numpy(float)
        ax = axes[0, j]
        ax.errorbar(p - 0.003, sm, yerr=np.vstack([sm - sl, sh - sm]), fmt="o", capsize=3, label="symmetric feasible domain")
        ax.errorbar(p + 0.003, am, yerr=np.vstack([am - al, ah - am]), fmt="s", capsize=3, label="asymmetric feasible domain")
        ax.axhline(1.0, linestyle="--", linewidth=1.0, label="2D Ising reference nu=1")
        ax.set_ylabel("Correlation-length exponent, nu")
        ax.set_title("Binder ratio-of-averages" if ch == "binder_roa" else "Correlation-length ratio xi/L")
        ax.grid(alpha=0.22)
        ax.legend(fontsize=7)

        ax2 = axes[1, j]
        ax2.plot(p, g.symmetric_boundary_fraction.to_numpy(float), "o-", label="symmetric")
        ax2.plot(p, g.asymmetric_boundary_fraction.to_numpy(float), "s--", label="asymmetric")
        ax2.axhline(max_boundary, linestyle=":", linewidth=1.0, label=f"locked gate={max_boundary:.2f}")
        ax2.set_ylabel("Boundary-hit fraction")
        ax2.set_xlabel("Site occupation probability, p")
        ax2.grid(alpha=0.22)
        ax2.legend(fontsize=7)

    for lab, ax in zip(labels, axes.ravel()):
        ax.text(0.01, 0.99, lab, transform=ax.transAxes, ha="left", va="top", fontweight="bold")
    fig.suptitle("Sensitivity to removing symmetry of the feasible nu domain", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fdir = out / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    pdf = fdir / "FIG_S1_ASYMMETRIC_NU_SENSITIVITY.pdf"
    png = fdir / "FIG_S1_ASYMMETRIC_NU_SENSITIVITY.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=400)
    plt.close(fig)
    meta = {
        "figure": "FIG_S1_ASYMMETRIC_NU_SENSITIVITY",
        "method_version": METHOD_VERSION,
        "meaning": "Presentation-only comparison of the original symmetric feasible-domain replay and the asymmetric-domain sensitivity, using paired bootstrap indices and unchanged locked scientific inputs.",
        "source_table": "tables/ASYMMETRIC_NU_RESULTS.csv",
        "source_table_sha256": sha256_file(out / "tables" / "ASYMMETRIC_NU_RESULTS.csv"),
        "pdf_sha256": sha256_file(pdf),
        "png_sha256": sha256_file(png),
    }
    atomic_json(fdir / "FIG_S1_ASYMMETRIC_NU_SENSITIVITY.metadata.json", meta)
    return meta


def output_manifest(out: Path) -> pd.DataFrame:
    rows = []
    target = out / "OUTPUT_SHA256_MANIFEST.csv"
    for p in sorted(out.rglob("*")):
        if not p.is_file() or p.resolve() == target.resolve():
            continue
        rows.append({"path": p.relative_to(out).as_posix(), "sha256": sha256_file(p), "bytes": int(p.stat().st_size)})
    df = pd.DataFrame(rows)
    atomic_csv(target, df)
    return df


def report_markdown(
    *,
    precheck: dict[str, Any],
    replay: dict[str, Any] | None,
    invariance: dict[str, Any],
    shadow: pd.DataFrame | None,
) -> str:
    lines = [
        "# Strict asymmetric-nu sensitivity audit",
        "",
        f"Method version: `{METHOD_VERSION}`",
        "",
        "This package is post-processing only. It does not alter Monte Carlo trajectories, the locked Tc estimator, bootstrap seeds/indices, support definition, interpolation, Pb residue, thresholds, lattice sizes, or the original decision tables.",
        "",
        f"Precheck: **{precheck['status']}**",
        f"Original-file invariance: **{invariance['status']}**",
    ]
    if replay is not None:
        lines.append(f"Symmetric replay against original production summaries: **{replay['status']}**")
    if shadow is not None:
        lines += ["", "## Shadow sensitivity statuses", ""]
        for _, r in shadow.sort_values("p").iterrows():
            lines.append(f"- p={r.p:.2f}: **{r.shadow_status}** - {r.reason}")
    lines += [
        "",
        "## Interpretation lock",
        "",
        "The shadow statuses are sensitivity diagnostics only. They never replace or overwrite the original v3.2.1/v3.2.1.1 scientific decision layer. If the symmetric replay fails, asymmetric results are not interpretable.",
    ]
    return "\n".join(lines) + "\n"


def load_user_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"config file not found: {path}")
    cfg = json.loads(path.read_text(encoding="utf-8"))
    if not cfg.get("project_root") or not cfg.get("spec_lock"):
        raise ValueError("USER_CONFIG.json must define project_root and spec_lock")
    return cfg


def main() -> int:
    ap = argparse.ArgumentParser(description="Strict post-processing sensitivity audit for asymmetric feasible nu domains.")
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--precheck-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    ns = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    cfg = load_user_config(ns.config.resolve())
    project_root = Path(cfg["project_root"]).expanduser().resolve()
    spec_path = Path(cfg["spec_lock"]).expanduser().resolve()
    workers = max(1, int(cfg.get("workers", 4)))
    out = (script_dir / "OUTPUT_ASYMMETRIC_NU_AUDIT").resolve()

    if not project_root.is_dir():
        raise SystemExit(f"FAIL-CLOSED: project_root does not exist: {project_root}")
    if not spec_path.is_file():
        raise SystemExit(f"FAIL-CLOSED: spec_lock does not exist: {spec_path}")
    spec_sha = sha256_file(spec_path)
    if spec_sha != EXPECTED_SPEC_SHA256:
        raise SystemExit(f"FAIL-CLOSED: SPEC_LOCK.json SHA-256 mismatch. Expected {EXPECTED_SPEC_SHA256}, got {spec_sha}")

    package_root = spec_path.parent.parent.resolve()
    package_src = package_root / "src"
    install_original_modules(package_src)
    spec = _LOAD_SPEC(spec_path)
    if spec.get("spec_name") != EXPECTED_SPEC_NAME:
        raise SystemExit(f"FAIL-CLOSED: unexpected spec_name: {spec.get('spec_name')}")

    original_results = (project_root / Path(spec["write_subdir"])).resolve()
    if not original_results.is_dir():
        raise SystemExit(f"FAIL-CLOSED: original results directory from locked spec does not exist: {original_results}")
    if is_within(out, project_root) or is_within(out, package_root) or is_within(out, original_results):
        raise SystemExit("FAIL-CLOSED: extract this sensitivity package outside the original project/audit-package trees. Output is not allowed inside original scientific directories.")
    out.mkdir(parents=True, exist_ok=True)

    paths, input_hashes = _VALIDATE_INPUTS(project_root, spec)
    cases = _LOAD_CASES(paths, spec)
    missing_labels = sorted(set(CASE_LABELS) - set(cases))
    if missing_labels:
        raise SystemExit(f"FAIL-CLOSED: required cases missing: {missing_labels}")

    result_paths = required_result_paths(original_results)
    source_paths = critical_source_paths(package_root)
    protected_paths = [spec_path, *source_paths, *paths.values(), *result_paths]
    before = snapshot(protected_paths)
    atomic_json(out / "manifests" / "BEFORE_SNAPSHOT.json", before)

    central = _CENTRAL_ALL(cases, spec)
    domain_map, task_bundle = build_domain_map(central, spec)
    atomic_csv(out / "tables" / "ASYMMETRIC_DOMAIN_MAP.csv", domain_map)

    selftest = run_selftests(spec)
    atomic_json(out / "manifests" / "SELFTEST_RESULT.json", selftest)
    if not selftest["passed"]:
        raise SystemExit("FAIL-CLOSED: asymmetric fitter self-test failed")

    method_lock = {
        "method_version": METHOD_VERSION,
        "created_utc": utc_now(),
        "spec_name": spec["spec_name"],
        "spec_sha256": spec_sha,
        "project_root": str(project_root),
        "original_results": str(original_results),
        "original_package_root": str(package_root),
        "scientific_change_under_test": "Remove only the forced symmetry of the centrally locked feasible nu domain about nu=1; retain all contiguous feasible segments within the declared domain.",
        "unchanged": [
            "Monte Carlo trajectories and input realization curves",
            "RG-invariant Tc estimator and per-bootstrap Tc re-estimation",
            "bootstrap base_seed and bootstrap_index sequence",
            "locked support target points",
            "Bhattacharjee-Seno Pb residue",
            "primary interpolation",
            "x-window and size windows",
            "declared nu bounds",
            "coarse-point count and scalar refinement tolerance",
            "boundary threshold 2 percent of the applicable feasible segment",
            "decision gates from SPEC_LOCK.json",
            "original result and decision files",
        ],
        "declared_nu_bounds": spec["nu_fit"]["declared_bounds"],
        "primary_interpolation": spec["nu_fit"]["primary_interpolation"],
        "coarse_points": int(spec["nu_fit"]["coarse_points"]),
        "refine_xatol": float(spec["nu_fit"]["refine_xatol"]),
        "feasibility_scan_points": SCAN_POINTS,
        "boundary_fraction_of_segment": BOUNDARY_FRACTION_OF_SEGMENT,
        "max_primary_boundary_fraction": float(spec["decision"].get("max_primary_boundary_fraction", 0.10)),
        "min_primary_valid_fraction": float(spec["decision"].get("min_primary_valid_fraction", 0.95)),
        "n_bootstrap": int(spec["bootstrap"]["n_large"]),
        "base_seed": int(spec["base_seed"]),
        "executor_backend": _BACKEND_NAME(),
        "input_hashes_from_original_validator": input_hashes,
    }
    atomic_json(out / "manifests" / "METHOD_LOCK.json", method_lock)

    precheck = {
        "status": "PASS",
        "passed": True,
        "spec_sha256_match": True,
        "original_input_validation": "PASS",
        "selftest": selftest["status"],
        "domain_rows": int(len(domain_map)),
        "protected_file_count": int(len(before)),
        "output_is_outside_original_trees": True,
    }
    atomic_json(out / "PRECHECK_RESULT.json", precheck)

    source_signature = {
        "method_version": METHOD_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "spec_sha256": spec_sha,
        "protected_before": {k: v["sha256"] for k, v in before.items()},
        "task_bundle_hash": stable_hash(task_bundle),
        "n_bootstrap": int(spec["bootstrap"]["n_large"]),
        "base_seed": int(spec["base_seed"]),
    }
    run_signature = stable_hash(source_signature)
    run_manifest_path = out / "manifests" / "RUN_MANIFEST.json"
    if run_manifest_path.exists():
        old = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        if old.get("run_signature") != run_signature:
            raise SystemExit("FAIL-CLOSED: existing OUTPUT_ASYMMETRIC_NU_AUDIT belongs to a different input/spec/source state. Rename/remove that output directory before continuing.")
    else:
        atomic_json(run_manifest_path, {**source_signature, "run_signature": run_signature, "created_utc": utc_now()})

    if ns.precheck_only:
        after = snapshot(protected_paths)
        atomic_json(out / "manifests" / "AFTER_SNAPSHOT.json", after)
        inv = compare_snapshots(before, after)
        atomic_json(out / "manifests" / "INVARIANCE_AUDIT.json", inv)
        atomic_text(out / "FINAL_ASYMMETRIC_NU_AUDIT_REPORT.md", report_markdown(precheck=precheck, replay=None, invariance=inv, shadow=None))
        output_manifest(out)
        if not inv["passed"]:
            raise SystemExit("FAIL-CLOSED: original protected files changed during precheck")
        print("PRECHECK PASS. No original scientific file was modified.")
        return 0

    boot = run_bootstrap(
        project_root=project_root,
        spec_path=spec_path,
        package_src=package_src,
        spec=spec,
        task_bundle=task_bundle,
        out=out,
        workers=workers,
        force=bool(ns.force),
        signature=run_signature,
    )
    summary = summarize_draws(boot, domain_map, spec)
    atomic_csv(out / "tables" / "ASYMMETRIC_NU_RESULTS.csv", summary)

    presentation_cols = [
        "case_label", "p", "channel",
        "symmetric_nu_median", "symmetric_nu_ci_low", "symmetric_nu_ci_high",
        "asymmetric_nu_median", "asymmetric_nu_ci_low", "asymmetric_nu_ci_high",
        "symmetric_boundary_fraction", "asymmetric_boundary_fraction",
        "symmetric_valid_fraction", "asymmetric_valid_fraction",
        "asymmetric_identifiability_gate_pass",
    ]
    atomic_csv(out / "tables" / "TABLE_PRESENTATION_SUMMARY.csv", summary[presentation_cols].copy())

    replay_df, replay = replay_audit(boot, summary, original_results)
    atomic_csv(out / "tables" / "SYMMETRIC_REPLAY_AUDIT.csv", replay_df)
    atomic_json(out / "manifests" / "SYMMETRIC_REPLAY_RESULT.json", replay)

    shadow = shadow_decisions(summary, spec, replay_pass=bool(replay["passed"]))
    atomic_csv(out / "tables" / "SHADOW_SENSITIVITY_DECISIONS.csv", shadow)
    atomic_json(out / "SHADOW_SENSITIVITY_RESULT.json", {
        "method_version": METHOD_VERSION,
        "scope": "sensitivity-only shadow result; original decision tables remain immutable",
        "rows": shadow.to_dict(orient="records"),
    })

    if replay["passed"]:
        make_figure(summary, out, spec)

    after = snapshot(protected_paths)
    atomic_json(out / "manifests" / "AFTER_SNAPSHOT.json", after)
    inv = compare_snapshots(before, after)
    atomic_json(out / "manifests" / "INVARIANCE_AUDIT.json", inv)

    final_pass = bool(precheck["passed"] and replay["passed"] and inv["passed"])
    final_result = {
        "status": "PASS" if final_pass else "FAIL",
        "passed": final_pass,
        "method_version": METHOD_VERSION,
        "precheck_pass": bool(precheck["passed"]),
        "symmetric_replay_pass": bool(replay["passed"]),
        "original_file_invariance_pass": bool(inv["passed"]),
        "shadow_decisions": shadow.to_dict(orient="records"),
        "interpretation_lock": "Do not modify the manuscript from these outputs unless symmetric replay and original-file invariance both PASS. Shadow decisions never overwrite original decisions automatically.",
    }
    atomic_json(out / "FINAL_RESULT.json", final_result)
    atomic_text(out / "FINAL_ASYMMETRIC_NU_AUDIT_REPORT.md", report_markdown(precheck=precheck, replay=replay, invariance=inv, shadow=shadow))
    output_manifest(out)

    if not inv["passed"]:
        raise SystemExit("FAIL-CLOSED: protected original files changed")
    if not replay["passed"]:
        raise SystemExit("FAIL-CLOSED: symmetric replay does not reproduce original production summaries; asymmetric results are not interpretable")
    print("FULL STRICT AUDIT PASS. See OUTPUT_ASYMMETRIC_NU_AUDIT/FINAL_RESULT.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
