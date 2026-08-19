# Scientific decisions locked for v3.2.1

1. **No new Monte Carlo is performed.** Accepted N60 realization curves are read-only inputs.
2. The legacy susceptibility-peak-shift Tc values obtained with fixed `nu=1` are **provenance and broad branch locators only**. Their confidence intervals are not optimizer bounds for the new nu inference.
3. The primary non-circular Tc anchor is obtained from adjacent-size crossings of two RG-invariant channels: Binder ratio-of-averages and `xi/L`. The largest two adjacent-size crossing estimates per channel are combined, then the two channel estimates are combined by their median.
4. The branch locator is deliberately broad (`±0.10`) and is audited at `±0.08` and `±0.12`; stability to this locator is reported.
5. `nu` is fit **after** the crossing Tc is determined, using a one-dimensional deterministic scalar search. There is no free 2-D `(Tc,nu)` optimizer valley in the primary inference.
6. `nu` feasibility is made symmetric around `nu=1` within the locked-support region. Candidates requiring extrapolation are invalid; locked points are never silently dropped.
7. Primary critical-scaling channels are Binder ratio-of-averages and `xi/L`. Free `q` fits are not inferential outputs.
8. `beta/nu` and `gamma/nu` are consistency diagnostics evaluated at the independently estimated crossing Tc; they are never used to tune Tc.
9. Quenched bootstrap resamples whole disorder-realization curves and preserves cross-observable dependence within a realization. The reported 95% interval is a **quenched-realization percentile-bootstrap interval**, not total uncertainty.
10. Every bootstrap draw re-estimates crossing Tc before estimating nu. Thus Tc uncertainty/covariance is propagated inside the quenched resampling path.
11. A paired `L_min` diagnostic compares `L=40..120` with `L=60..120` using the same bootstrap draw. Statistically significant motion toward `nu=1` vetoes a claim of asymptotic exponent change.
12. With only five sizes, no free correction exponent `omega` and no multi-parameter logarithmic correction law are fit. Finite-size/logarithmic corrections remain a competing explanation.
13. A separate robustness envelope combines bootstrap CI with conservative finite-size/Tc/interpolation/x-window sensitivities. **The robustness envelope is not a confidence interval.**
14. Exact pristine honeycomb-Ising theory is an external calibration gate only. It is not used to tune diluted cases. Calibration requires the exact pristine Tc, `nu=1`, `beta/nu=1/8`, and `gamma/nu=7/4` to lie within their corresponding correction-aware robustness envelopes.
15. Support imbalance is measured and reported. The published Bhattacharjee-Seno all-point residue is not post-hoc reweighted to obtain a preferred answer.
16. `EVIDENCE_AGAINST_NU1`, if ever emitted, is framework-specific evidence only. **A new universality-class claim is prohibited.**
17. No DFT exchange/anisotropy injection and no Kelvin material prediction is allowed from these accepted Ising trajectories.
18. Bootstrap checkpoints are atomic, signed, SHA-256 validated, resumable, and on Windows are placed in a short persistent workspace. A locked free-space reserve is checked before bootstrap begins.
