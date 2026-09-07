# SSLIP: review of the proposed reorganization

7 September 2026

Review of Philipp's [PR #2, “Implemented rotation, straightened out architecture”](https://github.com/Tijmenvermeij/SSLIP/pull/2), at commit [`222a9de`](https://github.com/PhilKro/SSLIP/commit/222a9de2448eb7e152a422d1b3bf0d335fd2bb06).

## Overall assessment

The proposed separation of preprocessing, solving, and plotting is worth adopting. It would make individual steps easier to test and allow figures to be regenerated without rerunning identification. The recommendation is to integrate this in small steps, because the PR also changes the public interface, output conventions, and interpretation of activity maps.

The points below concern the reviewed PR revision. They do not all apply to the current compatibility branch, which has already adapted selected fixes and rotation support. No further reorganization has been applied as part of this assessment.

## Which parts are useful?

| Proposal | Assessment |
| --- | --- |
| Separate preprocessing | Useful: filtering, coarse-graining, and gradient calculation can be tested and reused independently. |
| Separate solvers | Extracting the single-slip calculation is useful. Renaming the existing standalone solvers provides less benefit and requires callers to migrate. |
| Separate plotting | Useful, especially for deformation, residual, and rotation figures that can be regenerated from saved results. |
| Accept either displacements or gradients | Useful as a separate feature, with explicit rules for grid correspondence and whether preprocessing has already been performed. |
| Nested configuration | Reasonable internally; it does not require replacing the existing public options structure immediately. |
| Source subfolders and initialization | Potentially helpful after the functional separation, once example paths and dependencies are self-contained. |

## Issues to resolve before integration

### 1. Interface and example consistency

The [new entry point](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/SSLIP.m#L1) changes several conventions together:

- Five inputs become four, with `U` and `V` replaced by `DeformationData` and options split into configuration sections.
- The second output becomes preprocessed data instead of options.
- Activities change from systems-by-pixels to pixels-by-systems.

The supplied examples still use the old five-input call, which fails with “Too many input arguments.” Selecting a subset also leaves inconsistent indexing: requesting original system 2 returns a one-system list while retaining `NoSs = 2`, which then fails in plotting.

For the first integration, retain `[ebsdID, opt] = SSLIP(ebsd, U, V, sSLocal, opt)`, the current activity orientation, and original system numbering. The internal separation can be adopted behind that interface. Any later public interface change should include the examples and a documented migration.

### 2. Gradient inputs need a clear grid contract

In [preprocessSSLIP](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/preprocessSSLIP.m#L49), a gradient shape different from the original EBSD grid is effectively accepted as already processed. A small runtime check accepted 2-by-3 gradient arrays on a 4-by-4 spatial grid without an error.

All gradient components need compatible dimensions, spatial ordering, and a matching output grid. The code should distinguish raw and already processed inputs explicitly. With gradient-only input, the current implementation also inserts zero displacement fields. Those displacements are unavailable; plotting them as zero would imply information that was not supplied.

### 3. Activity grouping changes the physical interpretation

The [slip-activity plot](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_SlipActivities.m#L130) groups projected tensors that are approximately equal **or opposite**, then adds their signed activities without accounting for the tensor sign.

A runtime check with grouping enabled and no orientation supplied used two opposite tensors and activities `+0.2` and `-0.2`. Their deformation contributions add to `Hxy = 0.4`, but the plotted grouped activity was zero. Grouping needs to account for the sign relative to the representative tensor and define what the displayed quantity means.

Grouping and reporting logic are also duplicated between the ordinary and dominant-activity plots, with different default similarity thresholds (`0.05` and `0.2`). These analysis choices should be reviewed separately from extracting plotting functions.

### 4. Dominant-activity ranking can miss localized slip

The [dominant-activity plot](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L94) ranks systems using the 85th percentile of finite absolute activities, including zeros. A field with activity `0.2` on 10% of pixels and zero elsewhere therefore has a ranking value of zero. With the default threshold, that system goes into “leftovers.”

This is a meaningful choice for localized slip bands and needs a separate decision about the intended ranking. Changing the percentile alone is not necessarily a general solution.

### 5. Paths and dependencies need cleanup

`initSSLIP` adds `src` recursively, including `src/depr`, but does not add the relocated helpers in `examples/utils`. The examples also retain their older path setup. In addition, [dominant-plot export](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L371) calls `exportScaledFigure`, which is absent from the PR and the inspected official MTEX source; it resolved to a personal helper on the review laptop.

The examples and figure exports should run from a clean checkout with documented dependencies, without relying on another SSLIP checkout or personal helpers already on the MATLAB path.

## Proposed integration sequence

1. **Adapt the preprocessing separation and small plotting functions from this PR.** Retain the current public interface, solver names, output shapes, and numerical behavior. Build on Philipp's functions without introducing an additional framework.
2. **Verify this structural change on MTEX 6.1 and 7.0.** Compare each version against its own baseline, including activities, residuals, solver flags, selected-system labels, and examples with rotation off and on. Preserve zero activity for valid pixels below the identification threshold; failed or invalid fits remain distinguishable.
3. **Review the additional features separately.** Direct gradient input, grouping, dominant-system ranking, and automatic stress alignment each need their own behavior definition and targeted checks.
4. **Move folders and revise public interfaces only where useful.** Update initialization, examples, and dependency documentation together. Preserve the current sample dataset and paper citations.

Keep `coneprog` as the combined solver; no `fmincon` implementation is needed.

## Existing integration and contribution credit

Philipp's single-slip amplitude correction and explicit rotation basis have already been adapted on [`codex/mtex7-compatibility`](https://github.com/Tijmenvermeij/SSLIP/tree/codex/mtex7-compatibility), reviewed here at `09c8aac`. The adaptation also corrects normalized rotation units and guards unsupported solver methods. Both supplied examples run with rotation disabled and enabled under MTEX 6.1.0 and official 7.0.0. This establishes compatibility for those checks; it does not imply identical solver outcomes across MTEX versions.

The feature commits retain Philipp's contribution through source references and `Co-authored-by` trailers. Future adaptations should preserve the same clear attribution, distinguishing his original implementation from integration corrections.

The additional architecture findings above combine source review with small runtime probes under MATLAB R2024b and official MTEX 7.0.0. They are not a complete validation of the original PR on both MTEX versions.
