# Scientific scope lock

## Primary result that must not be replaced
The prespecified symmetric v3.2.1/v3.2.1.1 analysis remains the primary scientific decision layer. This package is a later hardening layer only.

## New computation allowed by this package
One post-hoc shadow sensitivity is allowed: rebuild the corrected full-window energy-ESS cell diagnostic and prohibit failed cells only when they would serve as cubic4 interpolation base/source points in the primary nu-collapse calculation.

Everything else stays locked: raw Monte Carlo tables, per-draw Tc, bootstrap indices, common target support, symmetric feasible nu bounds, lattice-size windows, primary cubic4 interpolation, Pb residue apart from the explicit source mask, and the original boundary/validity thresholds.

## Explicitly prohibited
- No new Monte Carlo.
- No regeneration or replacement of accepted N60 trajectories.
- No change of the 0.10 boundary-hit or 0.95 valid-fit gates.
- No widening/narrowing of nu bounds after observing the ESS result.
- No silent removal of locked target residuals.
- No fallback to linear/PCHIP when the masked cubic4 fit is invalid.
- No promotion of asymmetric or ESS sensitivities to the original prespecified analysis.
- No new-universality-class claim.
- No material-specific or Kelvin mapping.

## Full-window ESS terminology
The 100 / 90% full-window rule implemented here is a post-hoc sensitivity diagnostic. It is not the original publication near-Tc gate. The original near-Tc diagnostic and this full-window source-dependency sensitivity must remain clearly distinguished in any manuscript revision.

## Fail-closed manuscript rule
The package may mark manuscript support READY only when the diluted conclusion is stable and the package's explicit strong-closure conditions are met. Otherwise it emits HOLD_FOR_SCIENTIFIC_REVIEW and does not overwrite the manuscript.
