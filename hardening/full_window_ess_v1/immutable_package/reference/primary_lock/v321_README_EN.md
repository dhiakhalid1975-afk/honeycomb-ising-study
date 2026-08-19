# FGT Correction-Aware Critical Audit v3.2.1

A read-only post-processing audit for the accepted N60 site-diluted honeycomb-Ising data. The package does not rerun Monte Carlo and does not alter accepted trajectories.

Primary changes relative to v3.1.2:
- legacy fixed-nu susceptibility-shift Tc values are provenance/branch locators only;
- Tc for nu inference comes from Binder and xi/L crossings;
- each quenched bootstrap draw re-estimates Tc before nu;
- nu is a deterministic 1-D scalar fit, avoiding a free Tc-nu optimizer valley;
- paired Lmin drift is an inference veto for correction-dominated apparent exponents;
- exact p=1 Tc, nu, beta/nu and gamma/nu are external calibration gates;
- bootstrap percentile intervals and systematic robustness envelopes are explicitly separated;
- no free q inference, no free omega fit with five sizes, no DFT/Kelvin injection.

Windows run order:
`00_INSTALL_AND_VALIDATE_ENV.cmd` -> `01_SYNTHETIC_CHALLENGE_4_WORKERS.cmd` -> `02_BOOTSTRAP_CONVERGENCE_4_WORKERS.cmd` -> `03_RUN_OR_RESUME_REAL_AUDIT_4_WORKERS.cmd` -> `06_VALIDATE_RELEASE.cmd`.

Intermediate chunks are atomic, signed, SHA-256 checked and resumable. Figures are exported as PDF, SVG and 600-dpi PNG.
