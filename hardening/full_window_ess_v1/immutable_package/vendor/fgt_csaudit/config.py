from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .provenance import stable_json_hash


def package_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_spec_path() -> Path:
    return package_root() / "configs" / "SPEC_LOCK.json"


def load_spec(path: Path | None = None) -> dict[str, Any]:
    p = path or default_spec_path()
    spec = json.loads(p.read_text(encoding="utf-8"))
    if spec.get("spec_name") != "FGT_CORRECTION_AWARE_CRITICAL_AUDIT_v3.2.1":
        raise ValueError("unexpected or unlocked SPEC file")
    return spec


def spec_hash(spec: dict[str, Any]) -> str:
    return stable_json_hash(spec)
