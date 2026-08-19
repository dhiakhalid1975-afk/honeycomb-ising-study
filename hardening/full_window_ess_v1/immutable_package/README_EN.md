# HONEYCOMB_FINAL_HARDENING_v1_STRICT

Publication-support post-processing package for the quenched site-diluted honeycomb Ising study. It performs no new Monte Carlo and never overwrites the frozen v3.2.1/v3.2.1.1 primary decision layer.

The only new scientific computation is a post-hoc ESS source-dependency shadow sensitivity. Full-window corrected energy ESS is reconstructed from the original realization tables. Cells failing the locked diagnostic rule (ESS >= 100 in at least 90% of realizations) are prohibited only as cubic4 interpolation source/base cells. Locked targets, per-draw Tc, bootstrap indices, symmetric nu bounds, Pb definition, size windows, and decision thresholds remain unchanged. No interpolation fallback is allowed.

Run the numbered CMD files in order. The package is fail-closed and emits manuscript-ready tables/text only from completed outputs. It labels the full-window ESS rule as post-hoc and the completed asymmetric-domain analysis as a later rule-locked shadow sensitivity, not an original prespecification.
