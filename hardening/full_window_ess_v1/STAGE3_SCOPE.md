# Stage 3 - Full-window ESS hardening

Status: later post-hoc ESS source-dependency stress-test layer; not primary inference.
Immutable scientific package verified SHA256 PASS (72/72).
Final hardening outputs verified SHA256 PASS (20/20).
Excluded: checkpoints, caches, local machine-path files, USER_CONFIG, and embedded review ZIP.
The embedded review ZIP was excluded because the outer package manifest reports hash and size mismatch; this does not affect the verified 72-file immutable scientific layer or the verified 20-file output layer.
This layer must not overwrite the canonical v3.2.1.1 primary inference.
