NU_ASYMMETRIC_SENSITIVITY_v1_STRICT
===================================

PURPOSE
This package performs one controlled post-processing sensitivity test: it removes only the forced symmetry of the centrally locked feasible nu domain about nu=1. It uses the same locked v3.2.1 data, Tc estimator, bootstrap seed/index sequence, support, Pb residue, interpolation, size windows, declared nu bounds, and decision gates.

IT DOES NOT
- run new Monte Carlo;
- change any original project/audit file;
- overwrite TABLE_FINAL_DECISIONS_v321.csv;
- automatically change the manuscript conclusion.

STRICT SAFETY DESIGN
1. SPEC_LOCK.json must match the exact SHA-256 used by the final figure provenance:
   98f5b9bee3b4f39495d70c88ed0b922f102e12abb12df0324c0975bffb3aace2
2. Output must be outside the original project and original audit-package trees.
3. Critical source, input, and original result files are SHA-256 snapshotted before and after.
4. The original symmetric method is replayed on the same bootstrap indices.
5. Asymmetric results are interpretable only if that replay reproduces the original Tc/nu summaries within 5e-10 and original-file invariance passes.
6. The asymmetric feasible topology is locked centrally using the same 181-point feasibility scan used by the original symmetric-bound construction. Disconnected feasible segments are not bridged.

FILES TO RUN
- USER_CONFIG.json                 : only project_root/spec_lock/workers.
- RUN_01_PRECHECK_ONLY.cmd         : mandatory first step.
- RUN_02_FULL_STRICT_AUDIT.cmd     : full paired 2000-draw sensitivity run.
- RUN_03_PACK_OUTPUT_FOR_REVIEW.cmd : packages the completed OUTPUT folder for review/upload.
- RUN_02_FULL_STRICT_AUDIT_FORCE_REBUILD.cmd : use only if an intentional clean recomputation of checkpoints is required.

OUTPUT
All new files are written only under:
OUTPUT_ASYMMETRIC_NU_AUDIT

Key outputs after a successful full run:
- FINAL_RESULT.json
- FINAL_ASYMMETRIC_NU_AUDIT_REPORT.md
- tables/ASYMMETRIC_NU_RESULTS.csv
- tables/ASYMMETRIC_DOMAIN_MAP.csv
- tables/SYMMETRIC_REPLAY_AUDIT.csv
- tables/SHADOW_SENSITIVITY_DECISIONS.csv
- figures/FIG_S1_ASYMMETRIC_NU_SENSITIVITY.pdf
- manifests/INVARIANCE_AUDIT.json
- manifests/SYMMETRIC_REPLAY_RESULT.json
- OUTPUT_SHA256_MANIFEST.csv

INTERPRETATION LOCK
A shadow status is a sensitivity diagnostic only. It never replaces the original v3.2.1/v3.2.1.1 decision layer automatically.
