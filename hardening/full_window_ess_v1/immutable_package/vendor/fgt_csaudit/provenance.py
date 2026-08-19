from __future__ import annotations

import errno
import hashlib
import json
import os
import platform
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def safe_fsync(file_obj: Any) -> None:
    """Flush a regular file to stable storage where the platform supports it.

    On some Windows launch/filesystem combinations ``os.fsync`` can raise EBADF or
    EINVAL even for a file that has already been flushed successfully.  Atomic
    promotion and SHA-256 checkpoint validation remain the integrity mechanism in
    this package, so those Windows-specific fsync failures are treated as a
    durability limitation rather than a scientific/runtime failure.  All other
    errors remain fail-closed.
    """
    try:
        os.fsync(file_obj.fileno())
    except OSError as exc:
        if os.name == "nt" and exc.errno in {errno.EBADF, errno.EINVAL}:
            return
        raise


def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_json_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return sha256_bytes(payload)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".fgt_", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            safe_fsync(fh)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_write_json(path: Path, obj: Any) -> None:
    atomic_write_text(path, json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def environment_record() -> dict[str, Any]:
    import numpy as np
    import pandas as pd
    import scipy
    import matplotlib

    return {
        "created_utc": utc_now(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }


def atomic_write_csv(path: Path, df: Any) -> None:
    """Atomically write a pandas DataFrame without importing pandas at module import."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".fgt_", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False)
        # Flush file contents before promotion.
        with open(tmp, "rb") as fh:
            safe_fsync(fh)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def hash_tree(root: Path, *, suffixes: tuple[str, ...] | None = None) -> dict[str, str]:
    """Return deterministic relative-path -> SHA-256 map for a source tree."""
    out: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if suffixes is not None and p.suffix.lower() not in suffixes:
            continue
        rel = p.relative_to(root).as_posix()
        if any(part in {"__pycache__", ".pytest_cache", ".git", ".mypy_cache"} or part.endswith(".egg-info") for part in p.parts):
            continue
        out[rel] = sha256_file(p)
    return out
