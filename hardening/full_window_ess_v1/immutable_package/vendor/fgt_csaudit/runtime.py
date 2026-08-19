from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

from .provenance import atomic_write_text, stable_json_hash


def _external_workspace_requested() -> bool:
    return bool(os.environ.get("FGT_CSAUDIT_WORK_ROOT")) or os.name == "nt"


def _default_windows_work_root() -> Path:
    """Return a user-writable, persistent and deliberately short Windows work root.

    We avoid the deep scientific project path for bootstrap checkpoints because Windows
    path-length limits can be reached by nested checkpoint namespaces plus temporary
    filenames.  USERPROFILE is user-writable and persistent across reboots, which keeps resume
    semantics without requiring administrator privileges or registry changes.
    """
    override = os.environ.get("FGT_CSAUDIT_WORK_ROOT")
    if override:
        return Path(override).expanduser()

    # USERPROFILE is substantially shorter than a deep Desktop/project path, is
    # user-writable without elevation, and persists across reboots.
    userprofile = os.environ.get("USERPROFILE")
    if userprofile:
        return Path(userprofile) / "FGTCSA"

    local = os.environ.get("LOCALAPPDATA")
    if local:
        return Path(local) / "FGTCSA"

    # Fallbacks are only for unusual Windows environments.
    temp = os.environ.get("TEMP") or os.environ.get("TMP")
    if temp:
        return Path(temp) / "FGTCSA"

    return Path.cwd() / ".fgtcsa_work"


def checkpoint_workspace(project_root: Path, out_dir: Path, spec: dict[str, Any]) -> Path:
    """Resolve the stable checkpoint root for this project/software build.

    Scientific outputs remain under ``out_dir``.  Only resumable intermediate bootstrap
    chunks are redirected on Windows (or when FGT_CSAUDIT_WORK_ROOT is explicitly set).
    The workspace key is deterministic, so rerunning after interruption resumes exactly
    the compatible chunks; each chunk still carries and validates its run signature and
    SHA-256 hashes.
    """
    project_root = project_root.resolve()
    out_dir = out_dir.resolve()
    if not _external_workspace_requested():
        root = out_dir / "checkpoints"
    else:
        base = _default_windows_work_root()
        key = stable_json_hash({
            "project_root": str(project_root),
            "software_build_version": str(spec.get("software_build_version", "unknown")),
            "scientific_method_version": str(spec.get("scientific_method_version", "unknown")),
        })[:12]
        root = base / key / "cp"

    root.mkdir(parents=True, exist_ok=True)
    return root


def checkpoint_workspace_record(project_root: Path, out_dir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    root = checkpoint_workspace(project_root, out_dir, spec)
    # Mandatory atomic I/O probe.  A build must prove the runtime workspace can create,
    # fsync/flush, atomically promote, read, and delete a file before bootstrap begins.
    probe = root / ".io_probe"
    payload = "FGT_CSAUDIT_CHECKPOINT_IO_PROBE\n"
    atomic_write_text(probe, payload)
    got = probe.read_text(encoding="utf-8")
    if got != payload:
        raise RuntimeError("checkpoint workspace I/O probe read-back mismatch")
    probe.unlink(missing_ok=True)
    free_bytes = int(shutil.disk_usage(root).free)
    min_free = int(spec.get("runtime_safety", {}).get("minimum_checkpoint_free_bytes", 0))
    if min_free > 0 and free_bytes < min_free:
        raise RuntimeError(
            f"checkpoint workspace has only {free_bytes} free bytes; minimum locked safety reserve is {min_free}. "
            "Set FGT_CSAUDIT_WORK_ROOT to a drive with more free space and rerun; existing compatible chunks remain resumable."
        )

    return {
        "checkpoint_root": str(root),
        "checkpoint_free_bytes_at_probe": free_bytes,
        "minimum_checkpoint_free_bytes": min_free,
        "external_to_scientific_output": bool(root.resolve() != (out_dir / "checkpoints").resolve()),
        "override_env_present": bool(os.environ.get("FGT_CSAUDIT_WORK_ROOT")),
        "platform_os_name": os.name,
        "io_probe": "PASS",
        "resume_policy": "deterministic_workspace_plus_run_signature_and_sha256_chunk_validation",
        "scientific_outputs_remain_in": str(out_dir),
    }
