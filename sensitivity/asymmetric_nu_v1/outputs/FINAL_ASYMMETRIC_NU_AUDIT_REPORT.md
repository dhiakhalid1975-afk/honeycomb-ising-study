# Strict asymmetric-nu sensitivity audit

Method version: `NU_ASYMMETRIC_SENSITIVITY_v1_STRICT`

This package is post-processing only. It does not alter Monte Carlo trajectories, the locked Tc estimator, bootstrap seeds/indices, support definition, interpolation, Pb residue, thresholds, lattice sizes, or the original decision tables.

Precheck: **PASS**
Original-file invariance: **PASS**
Symmetric replay against original production summaries: **PASS**

## Shadow sensitivity statuses

- p=0.80: **ROBUST_NONIDENTIFIABLE_LIMITED_RANGE** - at least one primary channel still fails the locked boundary/validity identifiability gate under asymmetric domains
- p=0.85: **REQUIRES_FULL_DECISION_REVIEW** - both asymmetric primary gates pass and no significant paired drift veto remains; this sensitivity alone cannot replace the original decision layer
- p=0.90: **REQUIRES_FULL_DECISION_REVIEW** - both asymmetric primary gates pass and no significant paired drift veto remains; this sensitivity alone cannot replace the original decision layer
- p=1.00: **PRISTINE_DIAGNOSTIC_PASS** - both asymmetric primary gates pass and both nominal bootstrap intervals contain nu=1; calibration is not redefined here

## Interpretation lock

The shadow statuses are sensitivity diagnostics only. They never replace or overwrite the original v3.2.1/v3.2.1.1 scientific decision layer. If the symmetric replay fails, asymmetric results are not interpretable.
