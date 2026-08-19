from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

from .config import default_spec_path, package_root
from .pipeline_v321 import (
    initialize_run, ensure_central, run_synthetic, run_bootstrap_convergence,
    run_real_audit, status, rebuild_figures, validate_release, run_all,
)


def _path(s: str) -> Path:
    return Path(s).expanduser().resolve()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fgt-csaudit-v321",
        description="Correction-aware, Tc-unconditioned post-processing audit for accepted N60 data.",
    )
    p.add_argument("--project-root", required=True, type=_path, help="Root of FGT_Dilution_Study_Code_v2.4.0_PUBLICATION_STRICT")
    p.add_argument("--spec", type=_path, default=default_spec_path())
    p.add_argument("--workers", type=int, default=4)
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("validate-inputs")
    sub.add_parser("central")
    sub.add_parser("synthetic")
    q = sub.add_parser("bootstrap-convergence"); q.add_argument("--force", action="store_true")
    q = sub.add_parser("real"); q.add_argument("--force", action="store_true")
    sub.add_parser("status")
    sub.add_parser("figures")
    sub.add_parser("validate-release")
    q = sub.add_parser("all"); q.add_argument("--force", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")
    try:
        if args.command == "validate-inputs":
            spec, out, ctx = initialize_run(args.project_root, args.spec)
            result = {"status": "PASS", "output_dir": str(out), "run_signature": ctx["run_signature"], "core_audit": ctx["core_audit"], "runtime_workspace": ctx["runtime_workspace"]}
        elif args.command == "central":
            c = ensure_central(args.project_root, args.spec)
            result = {"status": "PASS", "output_dir": str(c["out"]), "central_fit_rows": int(len(c["central"]["fits"]))}
        elif args.command == "synthetic":
            result = run_synthetic(args.project_root, args.spec, workers=args.workers)
        elif args.command == "bootstrap-convergence":
            result = run_bootstrap_convergence(args.project_root, args.spec, workers=args.workers, force=args.force)
        elif args.command == "real":
            result = run_real_audit(args.project_root, args.spec, workers=args.workers, force=args.force)
        elif args.command == "status":
            result = status(args.project_root, args.spec)
        elif args.command == "figures":
            result = rebuild_figures(args.project_root, args.spec)
        elif args.command == "validate-release":
            result = validate_release(args.project_root, args.spec)
        elif args.command == "all":
            result = run_all(args.project_root, args.spec, workers=args.workers, force=args.force)
        else:
            raise AssertionError(args.command)
        print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
        return 0
    except KeyboardInterrupt:
        print("\nInterrupted. Completed atomic bootstrap chunks are preserved. Rerun the same command to resume.", file=sys.stderr)
        return 130
    except Exception as exc:
        msg = f"FAIL-CLOSED: {type(exc).__name__}: {exc}"
        print(msg, file=sys.stderr)
        tb = traceback.format_exc()
        print(tb, file=sys.stderr)
        try:
            log = package_root() / "RUNTIME_FAILURE_TRACEBACK.log"
            log.write_text(msg + "\n\n" + tb, encoding="utf-8")
            print(f"Diagnostic traceback saved to: {log}", file=sys.stderr)
        except Exception:
            pass
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
