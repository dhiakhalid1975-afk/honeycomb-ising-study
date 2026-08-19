from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Literal

import numpy as np
from scipy.interpolate import PchipInterpolator

Channel = Literal["abs_m", "chi_abs", "binder_roa", "xi_over_L"]
Interpolation = Literal["cubic4", "pchip", "linear"]


@dataclass(frozen=True)
class SupportLock:
    sizes: tuple[int, ...]
    target_indices: dict[int, np.ndarray]
    x_low: float
    x_high: float
    n_ordered_residuals: int
    reference_tc: float
    reference_nu: float
    x_window: tuple[float, float]
    edge_points: int

    def to_jsonable(self) -> dict:
        d = asdict(self)
        d["target_indices"] = {str(k): v.tolist() for k, v in self.target_indices.items()}
        return d


def channel_kind(channel: Channel) -> str:
    if channel == "abs_m":
        return "magnetization"
    if channel == "chi_abs":
        return "susceptibility"
    return "dimensionless"


def scale_xy(t: np.ndarray, y: np.ndarray, L: int, tc: float, nu: float, q: float, channel: Channel) -> tuple[np.ndarray, np.ndarray]:
    if tc <= 0 or nu <= 0:
        raise ValueError("tc and nu must be positive")
    x = (np.asarray(t, float) - tc) / tc * float(L) ** (1.0 / nu)
    yy = np.asarray(y, float)
    kind = channel_kind(channel)
    if kind == "magnetization":
        yy = yy * float(L) ** q
    elif kind == "susceptibility":
        yy = yy * float(L) ** (-q)
    return x, yy


def build_locked_support(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    nu: float,
    channel: Channel,
    q: float,
    x_window: tuple[float, float],
    edge_points: int = 2,
    minimum_points_per_size: int = 3,
) -> SupportLock:
    sizes = tuple(sorted(curves))
    if len(sizes) < 3:
        raise ValueError("at least three sizes are required")
    xref: dict[int, np.ndarray] = {}
    lows: list[float] = []
    highs: list[float] = []
    for L in sizes:
        t, y = curves[L]
        x, _ = scale_xy(t, y, L, tc, nu, q, channel)
        if np.any(np.diff(x) <= 0):
            raise ValueError(f"L={L}: temperatures/scaled x must be strictly increasing")
        if len(x) < 2 * edge_points + 2:
            raise ValueError(f"L={L}: too few temperatures for locked support")
        xref[L] = x
        lows.append(float(x[edge_points]))
        highs.append(float(x[-edge_points - 1]))
    lo = max(float(x_window[0]), max(lows))
    hi = min(float(x_window[1]), min(highs))
    if not lo < hi:
        raise ValueError("no common locked-support interval across sizes")

    target_indices: dict[int, np.ndarray] = {}
    total_target = 0
    for L in sizes:
        t, y = curves[L]
        x, _ = scale_xy(t, y, L, tc, nu, q, channel)
        idx = np.flatnonzero((x >= lo) & (x <= hi)).astype(int)
        if len(idx) < minimum_points_per_size:
            raise ValueError(f"L={L}: locked support has only {len(idx)} points")
        target_indices[L] = idx
        total_target += len(idx)
    n = (len(sizes) - 1) * total_target
    return SupportLock(sizes, target_indices, lo, hi, int(n), float(tc), float(nu), tuple(map(float, x_window)), int(edge_points))


def _interp_cubic4(xb: np.ndarray, yb: np.ndarray, xt: np.ndarray) -> np.ndarray:
    if len(xb) < 4:
        raise ValueError("cubic4 interpolation requires at least four base points")
    xt = np.asarray(xt, float)
    k = np.searchsorted(xb, xt, side="left")
    starts = np.clip(k - 2, 0, len(xb) - 4)
    inds = starts[:, None] + np.arange(4)[None, :]
    xs = xb[inds]
    ys = yb[inds]
    out = np.zeros(len(xt), dtype=float)
    for a in range(4):
        basis = np.ones(len(xt), dtype=float)
        for b in range(4):
            if b == a:
                continue
            denom = xs[:, a] - xs[:, b]
            if np.any(np.abs(denom) < 1e-15):
                raise ValueError("duplicate x values in cubic4 interpolation")
            basis *= (xt - xs[:, b]) / denom
        out += ys[:, a] * basis
    return out


def _interpolate(xb: np.ndarray, yb: np.ndarray, xt: np.ndarray, method: Interpolation) -> np.ndarray:
    if method == "cubic4":
        return _interp_cubic4(xb, yb, xt)
    if method == "pchip":
        return np.asarray(PchipInterpolator(xb, yb, extrapolate=False)(xt), float)
    if method == "linear":
        return np.interp(xt, xb, yb)
    raise ValueError(f"unsupported interpolation method: {method}")


def pb_residue(
    curves: dict[int, tuple[np.ndarray, np.ndarray]],
    *,
    tc: float,
    nu: float,
    q: float,
    channel: Channel,
    support: SupportLock,
    interpolation: Interpolation = "cubic4",
    invalid_penalty: float = 1e6,
) -> tuple[float, int, bool]:
    """Bhattacharjee-Seno P_b with residual exponent q_residual=1.

    The support is fixed before optimization. If a candidate parameter set would
    require extrapolation for any locked point, the candidate is invalid and is
    penalized rather than silently dropping that residual.
    """
    if tc <= 0 or nu <= 0 or q < 0 or not np.isfinite([tc, nu, q]).all():
        return float(invalid_penalty), support.n_ordered_residuals, False

    scaled: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for L in support.sizes:
        t, y = curves[L]
        x, yy = scale_xy(t, y, L, tc, nu, q, channel)
        if np.any(np.diff(x) <= 0):
            return float(invalid_penalty), support.n_ordered_residuals, False
        scaled[L] = (x, yy)

    residual_sum = 0.0
    n = 0
    for base_L in support.sizes:
        xb_all, yb_all = scaled[base_L]
        base_finite = np.isfinite(xb_all) & np.isfinite(yb_all)
        xb = xb_all[base_finite]
        yb = yb_all[base_finite]
        min_base = 4 if interpolation == "cubic4" else 2
        if len(xb) < min_base:
            return float(invalid_penalty), support.n_ordered_residuals, False
        for target_L in support.sizes:
            if target_L == base_L:
                continue
            xt_all, yt_all = scaled[target_L]
            idx = support.target_indices[target_L]
            xt = xt_all[idx]
            yt = yt_all[idx]
            # A non-finite locked target is a failed candidate; it is never dropped.
            if np.any(~np.isfinite(xt)) or np.any(~np.isfinite(yt)):
                return float(invalid_penalty), support.n_ordered_residuals, False
            # No extrapolation: locked points must remain bracketed for every candidate.
            if np.any(xt < xb[0]) or np.any(xt > xb[-1]):
                return float(invalid_penalty), support.n_ordered_residuals, False
            try:
                pred = _interpolate(xb, yb, xt, interpolation)
            except Exception:
                return float(invalid_penalty), support.n_ordered_residuals, False
            if np.any(~np.isfinite(pred)):
                return float(invalid_penalty), support.n_ordered_residuals, False
            residual_sum += float(np.sum(np.abs(yt - pred)))
            n += int(len(xt))
    if n != support.n_ordered_residuals:
        # This must never happen with locked support; treat as a software integrity failure.
        return float(invalid_penalty), n, False
    return residual_sum / float(n), n, True
