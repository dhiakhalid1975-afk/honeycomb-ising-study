from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Any, Callable


def backend_name() -> str:
    """Return the deterministic executor backend used by this package.

    Windows uses threads deliberately.  The scientific tasks are independent and
    deterministic, while ThreadPoolExecutor avoids Windows console/handle inheritance
    failures that can surface as OSError: [Errno 9] Bad file descriptor when a package
    is launched from CMD/Explorer.  POSIX systems retain process workers.

    The backend can be overridden only for software diagnostics with
    FGT_AUDIT_EXECUTOR=thread|process.  This variable is recorded by the caller where
    relevant; it never changes seeds, support, fits, or decision rules.
    """
    override = os.environ.get("FGT_AUDIT_EXECUTOR", "").strip().lower()
    if override in {"thread", "process"}:
        return override
    return "thread" if os.name == "nt" else "process"


def make_executor(
    workers: int,
    *,
    initializer: Callable[..., Any] | None = None,
    initargs: tuple[Any, ...] = (),
):
    n = max(1, int(workers))
    if backend_name() == "thread":
        return ThreadPoolExecutor(max_workers=n, initializer=initializer, initargs=initargs)
    return ProcessPoolExecutor(max_workers=n, initializer=initializer, initargs=initargs)
