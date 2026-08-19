from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

from .pb import SupportLock, pb_residue


@dataclass(frozen=True)
class NuFit:
    channel: str
    tc: float
    nu: float
    pb: float
    valid: bool
    boundary_hit: bool
    nu_bounds: tuple[float, float]
    nfev: int
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def restrict_curves(curves: dict[int, tuple[np.ndarray, np.ndarray]], sizes: list[int] | tuple[int, ...]) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    wanted = set(int(x) for x in sizes)
    out = {int(L): v for L, v in curves.items() if int(L) in wanted}
    if len(out) != len(wanted):
        missing = sorted(wanted - set(out))
        raise ValueError(f"missing sizes in curve restriction: {missing}")
    return out


def symmetric_feasible_nu_bounds(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    channel: str,
    support: SupportLock,
    spec: dict[str, Any],
    interpolation: str | None = None,
) -> tuple[float, float]:
    declared = tuple(float(x) for x in spec["nu_fit"]["declared_bounds"])
    method = str(interpolation or spec["nu_fit"]["primary_interpolation"])
    penalty = float(spec["nu_fit"]["invalid_penalty"])
    grid = np.linspace(declared[0], declared[1], 181)
    ok: list[float] = []
    for nu in grid:
        val, _, valid = pb_residue(
            curves, tc=float(tc), nu=float(nu), q=0.0, channel=channel,
            support=support, interpolation=method, invalid_penalty=penalty,
        )
        if valid and np.isfinite(val) and val < 0.1 * penalty:
            ok.append(float(nu))
    if not ok or min(ok) > 1.0 or max(ok) < 1.0:
        raise RuntimeError(f"{channel}: no feasible symmetric nu interval around 1")
    room = min(1.0 - min(ok), max(ok) - 1.0)
    if room <= 0.01:
        raise RuntimeError(f"{channel}: feasible nu room around 1 is too small: {room}")
    return 1.0 - room, 1.0 + room


def fit_nu_fixed_tc(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    channel: str,
    support: SupportLock,
    spec: dict[str, Any],
    nu_bounds: tuple[float, float],
    interpolation: str | None = None,
) -> NuFit:
    method = str(interpolation or spec["nu_fit"]["primary_interpolation"])
    penalty = float(spec["nu_fit"]["invalid_penalty"])
    lo, hi = map(float, nu_bounds)
    if not (lo < 1.0 < hi):
        raise ValueError("nu bounds must bracket 1")
    ncoarse = int(spec["nu_fit"]["coarse_points"])
    grid = np.linspace(lo, hi, ncoarse)

    def f(nu: float) -> float:
        return float(pb_residue(
            curves, tc=float(tc), nu=float(nu), q=0.0, channel=channel,
            support=support, interpolation=method, invalid_penalty=penalty,
        )[0])

    vals = np.asarray([f(float(x)) for x in grid], float)
    valid = np.isfinite(vals) & (vals < 0.1 * penalty)
    if not np.any(valid):
        return NuFit(channel, float(tc), float("nan"), float(penalty), False, True, (lo, hi), len(grid), "no valid coarse points")
    idx_valid = np.flatnonzero(valid)
    ibest = int(idx_valid[np.argmin(vals[idx_valid])])
    best_nu = float(grid[ibest]); best_pb = float(vals[ibest]); nfev = len(grid)

    left = max(0, ibest - 1); right = min(len(grid) - 1, ibest + 1)
    if left < right:
        a, b = float(grid[left]), float(grid[right])
        try:
            res = minimize_scalar(f, bounds=(a, b), method="bounded", options={"xatol": float(spec["nu_fit"]["refine_xatol"])})
            nfev += int(getattr(res, "nfev", 0))
            score = float(f(float(res.x)))
            nfev += 1
            if bool(res.success) and np.isfinite(score) and score < best_pb:
                best_nu, best_pb = float(res.x), score
        except Exception:
            pass
    width = hi - lo
    frac = 0.02
    boundary = (best_nu - lo) <= frac * width or (hi - best_nu) <= frac * width or ibest in {0, len(grid) - 1}
    return NuFit(channel, float(tc), best_nu, best_pb, True, bool(boundary), (lo, hi), int(nfev), "deterministic coarse grid + bounded scalar refinement")


def interpolate_observable_at_tc(curves: dict[int, tuple[np.ndarray, np.ndarray]], tc: float) -> tuple[np.ndarray, np.ndarray]:
    sizes: list[int] = []
    vals: list[float] = []
    for L in sorted(curves):
        t, y = curves[L]
        t = np.asarray(t, float); y = np.asarray(y, float)
        mask = np.isfinite(t) & np.isfinite(y)
        t, y = t[mask], y[mask]
        if len(t) < 3 or not (float(t.min()) <= float(tc) <= float(t.max())):
            continue
        order = np.argsort(t); t, y = t[order], y[order]
        val = float(PchipInterpolator(t, y, extrapolate=False)(float(tc)))
        if np.isfinite(val) and val > 0:
            sizes.append(int(L)); vals.append(val)
    return np.asarray(sizes, float), np.asarray(vals, float)


def exponent_ratio_at_tc(curves: dict[int, tuple[np.ndarray, np.ndarray]], tc: float, *, kind: str, sizes: list[int] | tuple[int, ...] | None = None) -> dict[str, Any]:
    c = curves if sizes is None else restrict_curves(curves, sizes)
    L, y = interpolate_observable_at_tc(c, float(tc))
    if len(L) < 3:
        return {"success": False, "ratio": float("nan"), "n_sizes": int(len(L))}
    x = np.log(L); z = np.log(y)
    slope, intercept = np.polyfit(x, z, 1)
    pred = intercept + slope * x
    resid = z - pred
    dof = max(1, len(L) - 2)
    s2 = float(np.sum(resid * resid) / dof)
    sx = float(np.sum((x - np.mean(x)) ** 2))
    slope_se = float(np.sqrt(s2 / sx)) if sx > 0 else float("nan")
    if kind == "magnetization":
        ratio, ratio_se = -float(slope), slope_se
    elif kind == "susceptibility":
        ratio, ratio_se = float(slope), slope_se
    else:
        raise ValueError(kind)
    return {
        "success": True,
        "ratio": ratio,
        "fit_se": ratio_se,
        "slope": float(slope),
        "intercept": float(intercept),
        "n_sizes": int(len(L)),
        "sizes": [int(x) for x in L],
    }
