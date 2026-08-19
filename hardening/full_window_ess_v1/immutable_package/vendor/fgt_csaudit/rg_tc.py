from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq


@dataclass(frozen=True)
class CrossingEstimate:
    channel: str
    estimate: float
    n_pairs_used: int
    selected_pairs: tuple[tuple[int, int, float], ...]
    all_pairs: tuple[tuple[int, int, tuple[float, ...]], ...]
    branch_center: float
    branch_half_width: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _roots_between_curves(
    t1: np.ndarray,
    y1: np.ndarray,
    t2: np.ndarray,
    y2: np.ndarray,
    *,
    grid_points: int,
) -> list[float]:
    t1 = np.asarray(t1, float); y1 = np.asarray(y1, float)
    t2 = np.asarray(t2, float); y2 = np.asarray(y2, float)
    m1 = np.isfinite(t1) & np.isfinite(y1)
    m2 = np.isfinite(t2) & np.isfinite(y2)
    t1, y1, t2, y2 = t1[m1], y1[m1], t2[m2], y2[m2]
    if len(t1) < 3 or len(t2) < 3:
        return []
    o1 = np.argsort(t1); o2 = np.argsort(t2)
    t1, y1, t2, y2 = t1[o1], y1[o1], t2[o2], y2[o2]
    lo = max(float(t1[0]), float(t2[0])); hi = min(float(t1[-1]), float(t2[-1]))
    if not lo < hi:
        return []
    f1 = PchipInterpolator(t1, y1, extrapolate=False)
    f2 = PchipInterpolator(t2, y2, extrapolate=False)
    grid = np.linspace(lo, hi, int(grid_points))
    d = np.asarray(f1(grid) - f2(grid), float)
    roots: list[float] = []
    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        da, db = d[i], d[i + 1]
        if not np.isfinite(da) or not np.isfinite(db):
            continue
        if da == 0.0:
            roots.append(a)
            continue
        if da * db > 0.0:
            continue
        try:
            r = float(brentq(lambda x: float(f1(x) - f2(x)), a, b))
        except Exception:
            continue
        if not roots or abs(r - roots[-1]) > 1e-8:
            roots.append(r)
    return roots


def estimate_crossing_tc(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    channel: str,
    branch_center: float,
    branch_half_width: float,
    largest_pairs_used: int = 2,
    adjacent_only: bool = True,
    grid_points: int = 4001,
) -> CrossingEstimate:
    sizes = sorted(int(x) for x in curves)
    pairs = list(zip(sizes[:-1], sizes[1:])) if adjacent_only else [
        (a, b) for i, a in enumerate(sizes) for b in sizes[i + 1:]
    ]
    all_rows: list[tuple[int, int, tuple[float, ...]]] = []
    selected: list[tuple[int, int, float]] = []
    lo_branch = float(branch_center) - float(branch_half_width)
    hi_branch = float(branch_center) + float(branch_half_width)
    for L1, L2 in pairs:
        roots = _roots_between_curves(*curves[L1], *curves[L2], grid_points=grid_points)
        all_rows.append((L1, L2, tuple(float(r) for r in roots)))
        in_branch = [r for r in roots if lo_branch <= r <= hi_branch]
        if in_branch:
            chosen = min(in_branch, key=lambda r: abs(r - float(branch_center)))
            selected.append((L1, L2, float(chosen)))
    selected.sort(key=lambda x: (x[1], x[0]), reverse=True)
    used = selected[:max(1, int(largest_pairs_used))]
    vals = np.asarray([x[2] for x in used], float)
    est = float(np.median(vals)) if len(vals) else float("nan")
    return CrossingEstimate(
        channel=str(channel), estimate=est, n_pairs_used=int(len(used)),
        selected_pairs=tuple(used), all_pairs=tuple(all_rows),
        branch_center=float(branch_center), branch_half_width=float(branch_half_width),
    )


def estimate_joint_tc(
    curves_by_channel: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]],
    *,
    branch_center: float,
    spec: dict[str, Any],
    sizes: tuple[int, ...] | list[int] | None = None,
    branch_half_width: float | None = None,
) -> dict[str, Any]:
    cfg = spec["tc_estimator"]
    channels = tuple(cfg["channels"])
    half = float(cfg["branch_half_width"] if branch_half_width is None else branch_half_width)
    out: dict[str, Any] = {}
    vals: list[float] = []
    for ch in channels:
        curves = curves_by_channel[ch]
        if sizes is not None:
            wanted = set(int(x) for x in sizes)
            curves = {L: v for L, v in curves.items() if int(L) in wanted}
        ce = estimate_crossing_tc(
            curves,
            channel=ch,
            branch_center=float(branch_center),
            branch_half_width=half,
            largest_pairs_used=int(cfg["largest_pairs_used"]),
            adjacent_only=bool(cfg["adjacent_pairs_only"]),
            grid_points=int(cfg["dense_grid_points"]),
        )
        out[ch] = ce.to_dict()
        if np.isfinite(ce.estimate):
            vals.append(float(ce.estimate))
    require_both = bool(cfg.get("require_both_channels_for_joint", True))
    if require_both and len(vals) != len(channels):
        joint = float("nan")
    elif not vals:
        joint = float("nan")
    elif str(cfg.get("joint_combine", "median")) == "median":
        joint = float(np.median(np.asarray(vals, float)))
    else:
        joint = float(np.mean(np.asarray(vals, float)))
    out["joint_tc"] = joint
    out["n_channels"] = int(len(vals))
    out["channel_spread"] = float(max(vals) - min(vals)) if len(vals) >= 2 else float("nan")
    out["branch_half_width"] = half
    return out
