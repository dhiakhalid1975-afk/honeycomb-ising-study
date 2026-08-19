from __future__ import annotations
import csv, hashlib, json, os, platform, sys, tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REFERENCE = PACKAGE_ROOT / "reference"
PRIMARY_LOCK = REFERENCE / "primary_lock"
ASYM_REF = REFERENCE / "asymmetric_completed"
CONFIG_DIR = PACKAGE_ROOT / "configs"

CASE_P = {"random_p080":0.80,"random_p085":0.85,"random_p090":0.90,"pristine_p100":1.00}
PRIMARY_CHANNELS = ("binder_roa","xi_over_L")
DILUTED = ("random_p080","random_p085","random_p090")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def sha256_tree(root: Path, suffix: str = '.py') -> str:
    h=hashlib.sha256()
    files=sorted(p for p in root.rglob(f'*{suffix}') if p.is_file())
    for p in files:
        rel=p.relative_to(root).as_posix().encode()
        h.update(len(rel).to_bytes(8,'big')); h.update(rel)
        data=p.read_bytes(); h.update(len(data).to_bytes(8,'big')); h.update(data)
    return h.hexdigest()


def stable_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=str(path.parent)); os.close(fd)
    Path(tmp).write_text(text,encoding='utf-8',newline='\n')
    os.replace(tmp,path)


def atomic_json(path: Path, obj: Any) -> None:
    atomic_text(path,json.dumps(obj,indent=2,ensure_ascii=False)+'\n')


def atomic_csv(path: Path, df) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=str(path.parent)); os.close(fd)
    df.to_csv(tmp,index=False)
    os.replace(tmp,path)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding='utf-8'))


def load_config(path: Path | None = None) -> dict[str,Any]:
    p=(path or (PACKAGE_ROOT/'USER_CONFIG.json')).resolve()
    cfg=load_json(p)
    override=os.environ.get('FGT_PROJECT_ROOT','').strip()
    if override: cfg['project_root']=override
    return cfg


def project_paths(project_root: Path) -> dict[str,Path]:
    base=project_root/'results'/'publication_strict_phase3'
    audit=base/'final_csaudit_v3_2_1_correction_aware'
    return {
      'fine':base/'tables'/'fine_realizations_n60.csv',
      'reference':base/'tables'/'reference_realizations.csv',
      'quality_gate':base/'quality_gate.json',
      'adaptive':base/'manifests'/'adaptive_final_decision.json',
      'tc_table':base/'publication_export'/'TABLE_TC_VS_P.csv',
      'locked_support':audit/'tables'/'LOCKED_SUPPORT_AND_SYMMETRIC_NU_BOUNDS.json',
      'spec_used':audit/'manifests'/'SPEC_LOCK_USED.json',
      'final_audit':audit/'FINAL_AUDIT_RESULT_v321.json',
    }



def verify_immutable_package_manifest() -> dict[str, Any]:
    """Verify the frozen code/reference payload while allowing local path/config edits."""
    manifest = PACKAGE_ROOT / "IMMUTABLE_PACKAGE_SHA256_MANIFEST.csv"
    if not manifest.exists():
        raise RuntimeError(f"immutable package manifest missing: {manifest}")
    import pandas as pd
    df = pd.read_csv(manifest)
    required = {"path", "sha256", "bytes"}
    if not required.issubset(df.columns):
        raise RuntimeError("immutable package manifest schema mismatch")
    failures = []
    for r in df.itertuples(index=False):
        q = PACKAGE_ROOT / str(r.path)
        if not q.exists():
            failures.append(f"missing:{r.path}")
            continue
        if int(q.stat().st_size) != int(r.bytes):
            failures.append(f"size:{r.path}")
            continue
        if sha256_file(q) != str(r.sha256):
            failures.append(f"sha256:{r.path}")
    if failures:
        raise RuntimeError("immutable package integrity failure: " + ", ".join(failures[:20]))
    return {"manifest": str(manifest), "files_verified": int(len(df)), "status": "PASS"}

def environment_record() -> dict[str,Any]:
    import numpy, pandas, scipy
    return {
      'created_utc':utc_now(),'python':sys.version,'platform':platform.platform(),
      'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__,
    }
