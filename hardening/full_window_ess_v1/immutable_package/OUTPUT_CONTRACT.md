# Output contract

A successful publication run produces under `OUTPUT_FINAL_HARDENING/`:

## Provenance / gates
- `manifests/PRECHECK_RESULT.json`
- `manifests/BASELINE_REPLAY_RESULT.json`
- `manifests/ESS_AUDIT_RESULT.json`
- `OUTPUT_SHA256_MANIFEST.csv`

## ESS audit and inference dependency
- `tables/ESS_CELL_AUDIT_REBUILT.csv`
- `tables/ESS_CELL_AUDIT_SUMMARY.csv`
- `tables/LOCKED_TARGET_ESS_AUDIT.csv`
- `tables/ESS_MASKED_DRAW_LEVEL.csv`
- `tables/INTERPOLATION_DEPENDENCY_MAP.csv`
- `tables/ESS_MASKED_NU_RESULTS.csv`
- `FINAL_ESS_DEPENDENCY_RESULT.json`
- `FINAL_ESS_DEPENDENCY_REPORT.md`

## Manuscript support (only after the scientific guard is evaluated)
- `manuscript_support/MANUSCRIPT_INTEGRATION_GUARD.json`
- `manuscript_support/MANUSCRIPT_HARDENING_PLAN.md`
- `manuscript_support/SUPPLEMENTARY_TABLE_S4_ASYMMETRIC.csv`
- `manuscript_support/SUPPLEMENTARY_FIGURE_S3_ASYMMETRIC.pdf/.png`
- `manuscript_support/SUPPLEMENTARY_TABLE_S5_ESS_DEPENDENCY.csv`

## Review archive
`RUN_06_PACK_OUTPUT_FOR_REVIEW.cmd` builds `HONEYCOMB_FINAL_HARDENING_REVIEW_OUTPUT.zip` from non-checkpoint outputs.

Checkpoint files are resumability artifacts and are intentionally excluded from the review archive.
