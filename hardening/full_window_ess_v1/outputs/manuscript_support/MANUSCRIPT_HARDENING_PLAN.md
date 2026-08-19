# Manuscript hardening insertion plan

Status: **READY FOR CONTROLLED MANUSCRIPT INTEGRATION**

This file does not overwrite the manuscript. The prespecified symmetric v3.2.1/v3.2.1.1 analysis remains primary.

## 1. Abstract
Replace the unqualified diluted-boundary sentence with wording equivalent to:

> Under the prespecified symmetric primary analysis, at every p < 1 at least one primary nu channel exceeds the boundary-hit threshold, and all six diluted upper bootstrap endpoints coincide with the symmetric feasible bound. An additional rule-locked asymmetric-domain sensitivity leaves the central estimates essentially unchanged; although boundary fractions fall below the diagnostic threshold at p = 0.85 and 0.90, Binder and xi/L never jointly exclude nu = 1.

Do **not** describe the asymmetric test as prespecified in the original study. The completed audit is a later rule-locked shadow sensitivity.

## 2. End of Sec. 4.3 — asymmetric method paragraph
> After completion of the primary symmetric-domain analysis, an additional rule-locked post-processing sensitivity analysis was specified before execution. This analysis removed only the forced symmetry of the centrally locked feasible nu domain about nu=1, while leaving the Monte Carlo data, Tc estimator, bootstrap draws, common-support definition, interpolation rule, collapse residue, lattice-size windows, declared search range 0.55<=nu<=1.45, and decision thresholds unchanged. The resulting outputs were treated as shadow sensitivity diagnostics and were not allowed to overwrite the primary v3.2.1/v3.2.1.1 decision layer.

## 3. Sec. 6.4 — asymmetric result
Use the generated Supplementary Table S4. The correct count is **5 of 6** diluted asymmetric upper endpoints at 1.45; the exception is p=0.90 Binder with upper endpoint 1.273515. Do not translate REQUIRES_FULL_DECISION_REVIEW into an identified exponent.

## 4. Sec. 6.6 — ESS wording
The current near-Tc numbers should be explicitly scoped as **the locked near-Tc probe temperature at L=120**. Do not present them as minima over the full simulated window.

Add the following only as a post-hoc sensitivity result, not an original gate:

> The post-hoc ESS source-dependency sensitivity did not overturn the conservative diluted-exponent conclusion. Failed full-window energy-ESS cells were prohibited only as cubic-interpolation source points, while the locked target support, per-draw Tc values, bootstrap indices, symmetric nu bounds, interpolation rule, and decision thresholds were unchanged. No diluted case produced joint Binder and xi/L evidence against nu=1 under this source mask.

### Diluted ESS source-mask result to report
Use Supplementary Table S5 for exact values. In the completed run, p=0.80 is numerically unchanged because no failed full-window ESS cell is used by the original cubic4 stencils. At p=0.85, the masked analysis remains limited by the locked boundary criterion. At p=0.90, both masked boundary gates pass, while both masked 95% intervals contain nu=1; therefore the source mask does not create joint evidence against nu=1.

The complete audit also retains a pristine post-hoc stress diagnostic. Under the later full-window mask, one pristine locked target cell fails that post-hoc cell rule and the masked pristine collapse becomes invalid. Do not use this later stress mask to redefine or reject the original pristine calibration: the full-window rule was not the original near-Tc gate. Keep this distinction explicit in repository documentation.

## 5. Sec. 7.2 and Conclusions
Where the manuscript says “Under the prespecified criterion”, change this to **“Under the prespecified symmetric primary criterion”**. State that the asymmetric-domain and ESS source-dependency checks are additional sensitivities and do not replace the primary decision layer.

## 6. Supplement numbering
- S5. Asymmetric-domain sensitivity
- Supplementary Table S4: generated `SUPPLEMENTARY_TABLE_S4_ASYMMETRIC.csv`
- Supplementary Figure S3: generated copy of the completed asymmetric figure
- S6. ESS-inference dependency sensitivity
- Supplementary Table S5: generated `SUPPLEMENTARY_TABLE_S5_ESS_DEPENDENCY.csv` (diluted cases only)
- Repository/reviewer decision summary: `FINAL_HARDENING_DECISION_SUMMARY.csv`

## 7. Submission-only items still required
Authors/affiliations/corresponding author and the public repository DOI or persistent URL remain submission metadata and are not invented by this package.
