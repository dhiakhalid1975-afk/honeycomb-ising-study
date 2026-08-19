# Package validation report

Status: **PASS — package construction, exact baseline replay, and the complete 2000-draw publication hardening sensitivity were executed successfully against the frozen project inputs recovered from the uploaded project archive.**

## Exact input lock
All strict input SHA-256 checks passed, including the accepted diluted N60 table, pristine reference table, locked support, SPEC_LOCK_USED, completed asymmetric outputs, and key vendored v3.2.1 source files. Exact hashes are recorded in `EXACT_PUBLICATION_RUN_PROVENANCE.json`.

## Baseline reproduction
- Local deterministic replay: 24/24 sentinel draw-channel checks PASS (bootstrap indices 0, 999, 1999; four cases; two primary channels).
- Maximum |delta nu| = 6.7724e-14.
- Maximum |delta Pb| = 3.7470e-16.
- The embedded completed asymmetric package also carries its full symmetric replay audit.

## Full-window corrected-energy ESS audit
The later full-window 100/90% rule is explicitly post-hoc and is **not** the original near-Tc publication gate. Rebuilt failed energy-cell counts are 13, 13, and 20 for p=0.80, 0.85, and 0.90; the pristine stress audit has 26 failed cells. No diluted locked target cell fails this post-hoc cell rule.

## Complete source-dependency sensitivity
The completed run contains 16,000 draw-channel rows (2000 draws x 4 cases x 2 primary channels). For the diluted cases:

- p=0.80: no failed full-window ESS cell is used by the original cubic4 stencils; masked results are numerically unchanged.
- p=0.85: failed-source dependency occurs in 0.7170 of Binder draws and 0.7105 of xi/L draws. The masked valid fraction is 0.999 in each channel, but the boundary fractions remain 0.4140 and 0.3505, so the locked identifiability gate still fails.
- p=0.90: failed-source dependency occurs in every draw. After source masking, Binder gives median nu=1.055825, 95% interval 0.874082-1.235000, boundary fraction 0.0505; xi/L gives median nu=1.082454, 95% interval 0.976114-1.234809, boundary fraction 0.0310. Valid fractions are 0.997 in both channels. Both intervals contain nu=1.

Therefore the diluted conservative conclusion is stable, and the package-level strong ESS-dependency closure flag is PASS **for the diluted cases only**: no diluted case produces joint Binder and xi/L evidence against nu=1 under the source mask.

## Pristine stress diagnostic
The later full-window source mask is deliberately not used to redefine the pristine primary calibration. One pristine locked target cell fails that post-hoc cell rule, and the source-masked pristine collapse is invalid. This is retained transparently in the full audit and is not reinterpreted as failure of the original near-Tc-calibrated primary analysis.

## Software validation
- Python syntax compilation: PASS.
- Unit/integrity tests after final manifest construction: **8/8 PASS**, including an independent four-point Lagrange interpolation cross-check against the masked Pb implementation.
- The earlier 100-draw smoke run remains software validation only and is superseded scientifically by the complete 2000-draw run above.

## Final packaging integrity
The review-archive output manifest excludes itself by design and was independently verified against every payload member. This packaging-only correction occurred after the completed scientific run; the scientific numerical modules and method lock are unchanged, as recorded in `EXACT_PUBLICATION_RUN_PROVENANCE.json`.
