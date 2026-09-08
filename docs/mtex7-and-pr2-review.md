# MTEX 7 validation and Philipp's PR #2

Date: 7 September 2026. Cleanup checkpoint: `f0ef743`.

## Decision

The initial compatibility update, through `09c8aac`, kept the existing public
interface and code organization. The regression checks and both examples run under MTEX 6.1.0 and the
[official MTEX 7.0.0 release](https://github.com/mtex-toolbox/mtex/releases/tag/mtex-7.0.0).
No broad MTEX API rewrite was needed for these paths.

Adopt Philipp's `singleSlipPerPixel` correction and explicit rotation basis in
the existing functions. Correct rotation output scaling and reject unsupported
solver methods. Keep the existing five-input call, activity array orientation,
and full slip-system list, avoiding the interface and indexing problems in the
original PR. That integration introduced no new library functions. Subsequent
reorganization is being applied in separate checkpoints, beginning with the
preprocessing separation described below. The review branch
`codex/philipp-pr2-review` preserves his original revision unchanged.

## Validation

Environment: MATLAB R2024b Update 2 on Apple Silicon, with Optimization,
Parallel Computing, and Image Processing Toolboxes. Final checks used separate
MATLAB sessions initialized with one MTEX version. Both sessions exited normally.

- MTEX 6.1.0: focused checks and both examples passed.
- Official MTEX 7.0.0 (release tag commit `4208a1ef6866676f3c21b15a50a55d7c426169a9`): same checks and examples passed.
- The laptop's MTEX 7 development checkout (`28b0d0f6b8ad947bff636099023b1280e13dc67e`) also passed the initial checks and examples.
- A real two-worker process pool passed the coneprog and per-pixel exit-flag check under official MTEX 7.0.0.
- Exported final activity plots from both official-release examples were visually inspected.

The tests cover positive and signed coneprog solves, fixed and relative residual
tolerances, the noise floor, skipped/invalid pixels, infeasible fits, single-slip
filtering and selection, residual plot edge cases, and all three trace layouts.
The example runner exercises filtering, coarse-graining, deformation plots,
activity plots, residual plots, and saving results.

### Numerical comparison: 6.1.0 versus official 7.0.0

| Measure | Ni superalloy | Virtual HCP |
| --- | ---: | ---: |
| Output pixels | 18,054 | 13,689 |
| Slip systems | 12 | 24 |
| Maximum input gradient difference | 4.94e-14 | 4.65e-16 |
| Maximum slip tensor difference | 0 | 5.00e-16 |
| Pixels with different solver success/failure outcomes | 3 | 2 |
| Maximum activity difference, commonly valid pixels | 1.50e-4 | 1.36e-4 |
| Maximum residual difference, commonly valid pixels | 5.29e-8 | 2.36e-8 |
| Maximum reconstructed gradient component difference | 2.30e-5 | 2.50e-5 |
| Maximum L1 objective difference | 3.25e-7 | 2.23e-7 |

These are runtime compatibility results, not a claim of identical activity
values or identical success masks. Direct coneprog reruns on the differing
pixels showed exit flag `1` versus `-7` after round-off-level input changes.
Flag `-7` means the search direction became small without satisfying the
requested convergence tolerances. All five differing fits were close to the
residual constraint boundary. Neither formulation nor solver tolerances were
changed to force agreement.

The Ni example produced 8,207/8,208 successful fits and 3/2 stalled fits for
6.1/7.0; HCP produced 3,309 successful fits and one stalled fit in each version,
at different pixels. The remaining pixels were skipped or invalid.

For method 1, `ebsdID.prop.solverExitFlag` now records the status of each solve.
`1` is success, other numeric values are coneprog exit conditions, and `NaN`
means no solve was attempted. Unsuccessful fits continue to return NaN activity
and residual. The direct solver's third output now contains these per-pixel
flags instead of the old zero-filled placeholder.

## Selected rotation integration

### Contribution provenance

Philipp's verified Git identity is
`PhilKro <56341570+PhilKro@users.noreply.github.com>`.
Both the single-slip correction and rotation support originate in his commit
[`d3ee1a5`](https://github.com/PhilKro/SSLIP/commit/d3ee1a54eb9f8248565e0063060e749f9f56b9ce)
and are present at the reviewed PR head
[`222a9de`](https://github.com/PhilKro/SSLIP/commit/222a9de2448eb7e152a422d1b3bf0d335fd2bb06).

The single-slip adaptation commit `235b21c` includes his `Co-authored-by`
trailer. The rotation adaptation `b64616d` uses the same trailer and names the
original commit. The source comments and README also credit him. These are adaptations
to the current functions, rather than unchanged cherry-picks of his larger
architecture commit. His original commits retain their authorship on the
review branch; no attribution in his history is rewritten.

The integration corrections are the normalization-to-radians conversion,
the method guard, and extraction of rotation only after the combined solve.
The progress count also includes pixels exactly at `minEeff`, matching the
existing solve condition when pure rotation is included with `minEeff = 0`.
The regression checks and per-pixel solver diagnostics were added during this
review. They are distinguished from Philipp's original feature contribution.

### Formulation and behavior

Philipp's basis is the small-angle displacement gradient
`Hrot = theta * [0 -1; 1 0]`, matching Eq. (11), Section 2.4, of
[Vermeij et al., Strain 61 (2025), e70000](https://doi.org/10.1111/str.70000).
Positive `theta` means `Hxy = -theta`, `Hyx = theta` for pure rotation; angles
are in radians. This is a linear small-angle correction, not a finite-rotation
iteration. The original slip formulation and L1 objective are described in
Eq. (5) of [Vermeij et al., Acta Materialia 243, 118502](https://doi.org/10.1016/j.actamat.2022.118502).
The supplied local copies were checked on pages 7 and 4-5, respectively.

The adopted code retains Philipp's optimization choices:

- Rotation is an additional basis in the combined `coneprog` fit, enabled by
  `opt.enableRotation = 1` (default 0).
- Positive slip constraints still allow both signs of rotation by introducing
  positive and negative rotation terms and subtracting their fitted amplitudes.
- All fitted coefficients, including rotation, retain unit weights in the L1
  objective. Without normalization this penalizes `sum(abs(slip)) + abs(theta)`.
  With normalization, the rotation penalty is `sqrt(2)*abs(theta)`, while the
  slip coefficients refer to the normalized in-plane tensors as before.
  The output conversion fixes units without changing the PR's objective.
- `minEeff` and the residual threshold/fraction keep their existing meanings.
  Pure rotation has zero effective strain and is skipped with a positive
  `minEeff`; set it to zero to fit these pixels. Solver flags distinguish
  skipped pixels (`NaN`) from solved pixels (`1`) and failed fits.

`SSLIP` returns the angle in `ebsdID.prop.rotationIDcor` as pixels-by-one;
`slipIDcor` stays systems-by-pixels, containing only the physical systems in
`opt.NoSs` order. Calling `SSLIPConeprogConstrMinAbs` directly with rotation
enabled returns the angle as its last activity row, following the PR's solver
convention. Its second and third outputs remain residuals and exit flags.
Methods 2 and 3 raise `SSLIP:RotationRequiresMethod1` before preprocessing if
rotation is requested. They cannot silently consume a physical slip channel.

The +SSLIP paper combines rotation correction with preselection, pairwise fits,
and additional acceptance criteria. Adopting this PR extension does not
implement that full workflow or establish unique identification with arbitrary
linearly dependent slip systems. No new rotation weighting or selection
algorithm was introduced.

### Validation after integration

The same existing functions and options work on both MTEX versions; no version
detection, compatibility wrapper, or alternate implementation was needed.

| Check, MATLAB R2024b | MTEX 6.1.0 | Official MTEX 7.0.0 |
| --- | --- | --- |
| Existing checks plus rotation cases | Passed | Passed |
| Full Ni superalloy example | Passed | Passed |
| Full virtual HCP example | Passed | Passed |
| Rotation-disabled example results versus prior same-version checkpoint | Exactly equal | Exactly equal |

The baseline comparison covers all saved activity, residual, solver-flag,
gradient, coordinate, and slip-tensor arrays, including NaNs. The cross-version
differences remain exactly those reported above: three changed solver outcomes
for Ni and two for HCP. Both final sessions exited normally and used
`-noFigureWindows -nosplash`.

Rotation checks cover both angle signs with signed or positive slip bounds,
normalized and unnormalized bases, pure slip, pure rotation, mixtures, a noisy
mixture, selected-system output and plotting, skipped/invalid pixels, and
unsupported-method errors. Residual assertions allow coneprog's actual default
constraint tolerance (1e-6); its fitting tolerances were not changed.

### Supplied examples

Both existing example scripts now explicitly expose
`IDoptions.enableRotation = 0`. Setting it to 1 fits the rotation present in
their data, adds a separate rotation plot in degrees, and adds `_rotation` to
the saved result name. The stored `rotationIDcor` remains in radians. Default
slip plots and zero activity at skipped low-strain pixels retain their existing
behavior; the separate comparison figures made during review are not part of
the example plotting code.

The HCP example's optional batch 2 previously claimed that `[3 14 24]` were its
active systems, although the generated steps use `[3 3 7 19]`. That subset now
comes directly from `unique(systems,'stable')`, giving `[3 7 19]`. The default
24-system batch is unchanged. The Ni example's grain-number filename fix was
already included in the earlier compatibility checkpoint.

Both examples were run with rotation off and on under MTEX 6.1.0 and official
7.0.0. The corrected HCP subset was also run with rotation enabled on both
versions. Rotation field dimensions, physical channel counts, distinct output
names, and plot exports passed. Representative Ni and HCP rotation plots were
visually checked. The saved default numeric results exactly match the previous
same-version checkpoints, including solver flags and NaNs. These checks used
`-noFigureWindows -nosplash`; no new library functions were introduced.

## Review of the original PR #2

Reviewed [PR #2](https://github.com/Tijmenvermeij/SSLIP/pull/2),
"Implemented rotation, straightened out architecture", at
`222a9de2448eb7e152a422d1b3bf0d335fd2bb06`.
Its merge base with current main is `1e2faa0` from 3 December 2025, so integration
must preserve subsequent main changes, including the newer example dataset.

### Confirmed issues in the original PR

1. **The supplied examples call an incompatible interface.** The new
   `SSLIP(ebsd, DeformationData, sSLocal, cfg)` accepts four inputs, while both
   examples still pass five. A direct legacy call fails with "Too many input
   arguments". The second output also changes from options to preprocessed data,
   and the result activity array is transposed. A migration or compatibility
   layer is needed. See [the new entry point](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/SSLIP.m#L1).

2. **Selecting a subset leaves inconsistent system indices.** With two input
   systems and `cfg.solver.NoSs = 2`, the returned slip-system list has length
   one, while `optOut.NoSs` remains 2. Passing the returned objects to the slip
   plotting function raises an indexing error. Distinguish original system
   labels from indices into the selected list. See [selection](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/SSLIP.m#L88)
   and [plotting](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_SlipActivities.m#L81).

3. **Rotation silently consumes a physical slip channel in methods 2 and 3.**
   Only the combined solver appends rotation terms, but the common result mapping
   always removes the last activity channel when `enableRotation` is true.
   A two-system test returned only one physical-slip column in methods 2 and 3.
   Guard unsupported combinations or implement rotation consistently across
   solvers. See [result mapping](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/SSLIP.m#L114).

4. **Normalized rotation output is not rescaled to its original units.**
   A pure small rotation of 0.1 with `normalizeInplane = 1` returned 0.141421.
   Normalizing the antisymmetric basis divides it by sqrt(2), so its fitted
   coefficient must be converted back before being interpreted as an angle.
   See [normalization](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/solveSSLIP_Combined.m#L45).

### Useful work and remaining design choices

- The relative residual threshold and minimum noise floor already match the
  retained local cleanup, including the single-slip variable-name correction.
- The single-slip selection fix is ported: preserve the original amplitudes
  while creating the sparse selected output. Regression tests reproduce the
  old zero-output bug and verify both possible winning systems.
- Separating preprocessing, solvers, and plotting could make the library easier
  to extend, but is not necessary for MTEX 7 compatibility. Introduce it with
  updated examples and an explicit decision about public interfaces and array
  orientation.
- Explicit rotation control is adopted instead of overloading `CRSS == 0`.
  The PR's L1 penalty and low-effective-strain cutoff are retained and specified
  above; this is the selected rotation extension, not the full +SSLIP workflow.
- Preserve the current sample dataset and citation information during a future
  merge. The PR removes the sample MAT file and removes the paper citation
  header from the combined solver.

## Reorganization checkpoint 1: preprocessing

The [architecture assessment](philipp-pr2-reorganization-review.md) recommends
separating structural changes from new analysis features. This first checkpoint
adapts Philipp's `preprocessSSLIP` from `d3ee1a5` into the existing `src/` folder.
`SSLIP` now calls it for grid mapping, displacement filtering, coarse-graining,
gradient calculation, and storage of the processed fields.

The helper retains the current gridify step so displacement values follow their
coordinates, including when the input EBSD points are unordered. It accepts
Philipp's displacement structure and preprocessing options, with defaults still
supplied by `SSLIP`. Direct gradient input remains a separate proposed feature.
Dimension validation now rejects a mismatch in either input dimension; the old
vector-valued condition could miss a mismatch in only one dimension.

The public five-input call, returned options, solver functions, plotting code,
activity dimensions, and original slip-system numbering are retained. The
existing convention that zero displacement components become NaN when filtering
is enabled is also preserved; this is separate from zero fitted activity below
`minEeff`. Source comments identify Philipp's contribution, and the adaptation
commit includes his `Co-authored-by` trailer.

Validation used MATLAB R2024b in separate sessions with figure windows disabled.

| Check against the pre-reorganization checkpoint | MTEX 6.1.0 | MTEX 7.0.0 |
| --- | --- | --- |
| Focused checks, including preprocessing | Passed | Passed |
| Ni example, rotation off | Exactly unchanged | Exactly unchanged |
| Ni example, rotation on | Exactly unchanged | Exactly unchanged |
| HCP example, rotation off | Exactly unchanged | Exactly unchanged |
| HCP example, rotation on | Exactly unchanged | Exactly unchanged |

The exact comparisons use each MTEX version's own baseline. They cover all
saved activity, residual, solver-flag, gradient, coordinate, and slip-tensor
arrays, plus rotation when enabled and the returned options. NaNs and zeros are
included. This does not change the small cross-version differences documented
above. Both examples also exported their activity figures and, when enabled,
rotation figures successfully.

New focused checks use known affine displacements on a shifted rectangular
grid, shuffled EBSD points, coarse-graining, zero/missing displacement masks,
and an input mismatch in only one dimension. These use square pixels; support
for unequal pixel spacings under MTEX 6.1 is not established by this checkpoint.

## Reorganization checkpoint 2: small plotting routines

This checkpoint adapts Philipp's `plotSSLIP_DeformationFields`,
`plotSSLIP_Residual`, and `plotSSLIP_Rotation` from `d3ee1a5` into the existing
`src/` folder. `SSLIP` calls the deformation plotter, `plotSSLIP` calls the
residual plotter, and both examples call the rotation plotter. They can also
be called directly using the saved `ebsdID` and `optOut` variables. The rotation
plotter returns a figure handle for export, or an empty handle when no rotation
field is present.

The adaptation preserves the figures from checkpoint 1 (`c4d1134`):

- Deformation plots retain their layout, linear/logarithmic effective-strain
  scale, gradient limits, colormaps, font/size settings, and export filenames.
- Residual plots retain the full finite residual range and the non-degenerate
  fallback for zero or missing data. When `residualScaleSame` is supplied,
  `plotSSLIP` explicitly passes the activity limits as `opt.caxisMinMax`.
  A direct call with this option also requires those limits. The explicit
  residual argument remains authoritative even if the EBSD object already
  contains a different residual field.
- Rotation plots retain the examples' degree units, blue-to-red colormap,
  symmetric finite-maximum limits, labels, and lowercase `_rotation.png` name.
  Stored angles remain in radians, including zero values for skipped pixels.

No activity grouping, ranking, percentile clipping, context overlays, folder
moves, or solver changes are included. These functions are adaptations of
Philipp's separation using the existing plotting behavior; his contribution
is identified in source comments and the commit's `Co-authored-by` trailer.

| Check against checkpoint 1 | MTEX 6.1.0 | MTEX 7.0.0 |
| --- | --- | --- |
| Focused checks, including standalone plotters and exports | Passed | Passed |
| Ni and HCP numerical results, rotation off and on | Exactly unchanged | Exactly unchanged |
| Ni deformation, residual, and rotation figures | Pixel-identical | Pixel-identical |
| HCP deformation, residual, and rotation figures | Pixel-identical | Pixel-identical |

Numerical comparisons include all saved arrays and returned options, with
zeros and NaNs included, against each version's own baseline. The 12 figure
comparisons use saved rotation-enabled example results. They compare rendered
PNG pixels as well as titles, color scales and limits, colormaps, axis positions,
patch values/geometry, and colorbar labels/ticks. Representative regenerated
Ni deformation and rotation figures were also inspected visually.

Focused checks additionally exercise custom deformation limits, logarithmic
shared residual scales, stale stored residuals, unchanged radians after plotting,
missing rotation fields, and the existing PNG/JPEG export names. Both MATLAB
sessions exited successfully and used `-noFigureWindows -nosplash`.

## Reorganization checkpoint 3: single-slip solver

This checkpoint adapts Philipp's `solveSSLIP_SingleSlip` from `d3ee1a5` into
`src/`. The method-3 branch in `SSLIP` now delegates to this function. It still
uses the existing `SSLIPConstr` least-squares fit for each selected system;
the existing method-1 and method-2 solver functions retain their names and code.

The helper takes the already selected systems in their supplied order and
returns systems-by-pixels activities and residuals. It does not apply `NoSs`
again. An optional third output returns the executed options, preserving the
public `SSLIP` behavior that reports `posConstr = 0` after the single-slip
override. Rotation remains unsupported for method 3, including a direct call
to the extracted helper.

The existing fitting conventions are retained:

- Acceptance uses a strict residual comparison, with the same absolute floor
  and optional relative threshold.
- `singleSlipPerPixel` retains the winning amplitude and resolves an exact tie
  using the first system in the supplied order.
- Rejected activities become zero while their per-system residuals remain
  available. A NaN-gradient test retains zero activity and NaN residual.
- Method 3 retains its residual-based acceptance and does not apply method 1's
  `minEeff` cutoff. This is pre-existing behavior, not a new change to the
  combined solver's handling of skipped pixels.

| Check against checkpoint 2 (`d074c15`) | MTEX 6.1.0 | MTEX 7.0.0 |
| --- | --- | --- |
| Focused checks, including the extracted single-slip path | Passed | Passed |
| Ni single-slip results, six option combinations | Exactly unchanged | Exactly unchanged |
| HCP single-slip results, six option combinations | Exactly unchanged | Exactly unchanged |
| Both full examples, rotation off and on | Exactly unchanged | Exactly unchanged |

The 24 single-slip comparisons cover full and reordered subsets, fixed and
relative thresholds, one-system-per-pixel selection, normalization, and the
positive-constraint override. They reuse the saved processed displacements
with filtering disabled and coarse-graining set to 1, on a fresh dummy grid
carrying the saved coordinates. The supplied local slip systems are unchanged.
Activities, residuals, coordinates, original system labels, and returned
options exactly match each MTEX version's own pre-extraction baseline.

The full method-1 examples also exactly retain their saved numerical arrays
and returned options, including rotations, zeros, NaNs, and solver flags.
Additional focused checks cover known signed amplitudes, exact ties, normalized
coefficients, matrix-shaped gradient input, and strict threshold rejection.
Both MATLAB sessions exited successfully with figure windows disabled.

Source comments and the adaptation commit's `Co-authored-by` trailer identify
Philipp's contribution. The conservative separation of preprocessing, the small
plotting routines, and single-slip solving is now in place. Direct gradient
input, stress alignment, activity grouping/ranking, and broader interface or
folder changes remain separate proposals.

## Feature checkpoint 4: direct gradient input

This checkpoint adapts Philipp's structured displacement/gradient input from
`d3ee1a5`, alongside the existing five-input call:

```matlab
[ebsdID,optOut] = SSLIP(ebsd,U,V,sSLocal,opt);
[ebsdID,optOut] = SSLIP(ebsd,deformationData,sSLocal,opt);
```

The second output remains options, the options structure remains flat, and
activities retain their systems-by-pixels shape and original system labels.
A structure containing `U,V` takes the existing displacement path. A structure
containing `Hxx,Hxy,Hyx,Hyy` supplies displacement gradients directly.

### Input contract

The first gradient-input implementation accepts fields that are ready for
fitting. It does not infer whether they are raw or processed from their shape.
All four components must match the size and point order of the supplied EBSD
data. Each component is attached before gridify, so unordered points retain
their correspondence with the coordinates. A processed gradient field therefore
requires its matching processed grid, not the original higher-resolution grid.

Gradient input defaults to `filterSize = 0` and `coarsegrain = 1`; other values
are rejected explicitly. Displacement defaults remain unchanged. Mixed input,
incomplete components, mismatched dimensions, complex values, and infinities
are rejected. Real single/double arrays are accepted and converted to double;
zeros retain their physical meaning and NaN denotes missing data.

The output starts with clean data properties on a copy of the supplied grid,
retaining its spatial/orientation metadata without stale fits. Gradient-only
results omit `U,V` rather than reconstructing or filling unavailable
displacements. Their deformation figure contains the five available fields.
Zero or missing gradient fields receive non-degenerate display limits; stored
values and solver acceptance rules are unchanged.

Both examples retain the displacement call and now include a commented example
for refitting their saved gradients, with explicit preprocessing settings and
a separate case name. No new library functions were introduced for this feature.

### Validation

| Check | MTEX 6.1.0 | MTEX 7.0.0 |
| --- | --- | --- |
| Focused checks, including structured and direct-gradient input | Passed | Passed |
| Existing Ni/HCP examples, rotation off and on, versus checkpoint 3 | Exactly unchanged | Exactly unchanged |
| Ni/HCP refits from their saved gradients, rotation off and on | Exactly match existing fits | Exactly match existing fits |
| Existing deformation/residual/rotation figures | Pixel-identical | Pixel-identical |

The eight direct-gradient refits match the existing activity, residual,
solver-flag, gradient, effective-strain, and optional rotation arrays exactly,
including zeros and NaNs, against each version's own results. Gradient-input
options correctly report no additional filtering or coarse-graining; the
original displacement-path returned options remain exactly unchanged.

Focused checks cover reordered points with spatially varying components,
physical zeros and NaNs, malformed inputs, rejection of extra preprocessing,
known rotation/slip coefficients, removal of stale rotation and solver flags,
zero-field plots, and equivalence of structured and five-input displacements.
The 12 existing figure comparisons retain their graphics properties and PNG
pixels. Additional checks cover single-precision conversion and all-missing
gradient plots. The new Ni gradient-only figures were exported for inspection.
All runs used `-noFigureWindows -nosplash` in separate MTEX sessions.

This implements a bounded version of Philipp's input feature. Filtering raw
gradient fields, automatic stress alignment, and activity grouping/ranking
remain separate proposals. Source references and the feature commit's
`Co-authored-by` trailer retain Philipp's contribution credit.

## Feature checkpoint 5: stress alignment for positive activity

This checkpoint adapts Philipp's stress-alignment block from
[`9c393a3`](https://github.com/PhilKro/SSLIP/commit/9c393a3b41951f546cb1c04a95cd1ef294053176).
For methods 1 and 2 with `posConstr = 1`, `opt.stress` determines the signs of
the slip directions. A single `stressTensor` is interpreted in specimen
coordinates; a single nonzero `vector3d` specifies uniaxial tension.
The physical contraction of stress with each deformation tensor supplies the
resolved-shear sign. Negative signs flip `b`; exactly zero signs leave the
system intact. There is no Schmid-factor or CRSS normalization of this sign
calculation, and stress magnitude does not weight the objective.

The block is integrated directly into `SSLIP`, with flat options and no new
library function or version-specific path. It validates the stress type and
rejects nonfinite values, lists of stresses, and zero tension directions.
Without stress, it retains Philipp's warning and uses the caller's directions.
Signed fits do not align systems. Method 3 also skips alignment because it
switches `posConstr` off; flipping directions beforehand would needlessly
change its signed output. Rotation remains signed under positive slip bounds.

### Returned systems and adapted examples

The fitting basis is now available as a third output:

```matlab
[ebsdID,optOut,sSLocal] = SSLIP(ebsd,U,V,sSLocal,opt);
```

This is the **full column list**, including unselected systems, after alignment.
Activity rows still follow `optOut.NoSs`, which refers to the original labels.
Returning the full list avoids the original PR's selected-list/index mismatch.
Both displacement and structured-input calls expose the same output.

Both examples explicitly provide their uniaxial tensile stress and save the
returned systems with the numerical results. They retain their original signed
fit by default. The HCP example still sets its slip signs before generating
synthetic displacements; the stress-alignment step applies when fitting.

The saved list now has column orientation. Consequently, its numeric tensor
array has size `3 x 3 x N` rather than `3 x 3 x 1 x N`. A first strict storage
comparison caught this expected change; comparing each physical `3 x 3` tensor
confirmed identical entries and system order. No row-orientation compatibility
adapter was added. Returned example options now explicitly include `stress`.

### Validation

The focused checks pass on MTEX 6.1.0 and official 7.0.0. They exercise a general
stress tensor, uniaxial vector/tensor equivalence, reordered subsets, unchanged
CRSS, normalized coefficients, rotation, and reconstruction from the returned
systems for methods 1 and 2. Additional checks cover zero resolved shear, zero
stress, signed fitting, method 3, malformed inputs, and the missing-stress warning.

Both complete Ni and HCP examples, with rotation off and on, retain their
activities, residuals, solver flags, gradients, coordinates, and optional
rotation arrays exactly against each MTEX version's checkpoint-4 baseline.
The physical tensors also match exactly after accounting for list orientation;
the only option change is the explicitly supplied stress.

Positive-constrained refits of both examples on MTEX 7, with rotation off and
on, also exactly match fits using their existing manually oriented systems
without automatic alignment. The returned physical tensors are identical.

A controlled 63-pixel example on MTEX 7 imposes negative xy shear and a signed
rotation field. With the original slip direction and positive activity bounds,
no pixel has a feasible fit. After stress alignment, all 63 pixels solve
successfully. The maximum slip-activity error is `1.51e-6`; the maximum rotation
error is `1.05e-6` radians. These are known-input checks of a small independent
basis, not a claim of unique slip identification in the full examples.

All MATLAB runs disable figure windows. MTEX 7's local startup reports a
help-search indexing error after setting its paths and preferences; the SSLIP
checks and example calculations themselves complete successfully.
Source references and the feature commit's `Co-authored-by` trailer credit
Philipp. Activity grouping/ranking and folder organization remain separate steps.

## Reorganization checkpoint 6: folders and initialization

This checkpoint adopts Philipp's folder layout and `initSSLIP` from
[`d3ee1a5`](https://github.com/PhilKro/SSLIP/commit/d3ee1a54eb9f8248565e0063060e749f9f56b9ce).
The ten relocated functions retain their exact bytes, including existing author
and third-party notices. Solver names and numerical implementations stay intact.

| Folder | Contents |
| --- | --- |
| `src` | Entry point, preprocessing, effective strain, filtering, and solvers |
| `src/plotting` | The four existing activity/deformation/residual/rotation plotters |
| `src/utils` | `coarsegrainDisp`, `dummyEBSDSimple`, `nanconv`, and `parfor_wait` |
| `examples/utils` | `generateSlipStepField` and `getTraceProjB` |

The root-level initializer explicitly adds the three library folders. This
adapts Philipp's recursive path setup so deprecated, test, and example folders
are not automatically exposed as library code. It does not initialize MTEX,
change the working directory, or save the user's MATLAB path.

Each example derives the repository root from its own filename, calls
`initSSLIP`, and adds `examples/utils`. Ni data lookup now also uses that root.
This resolves the PR's missing example-helper paths and removes dependence on
`pwd` for input lookup. Outputs continue to use the execution directory; the
README explains MATLAB's `run` behavior and the isolated example runner.

For existing analysis scripts, replace a path setup that only adds `src` with:

```matlab
addpath('/path/to/SSLIP');
initSSLIP;
```

No forwarding functions were left in the old locations. The existing examples
and checks initialize the new layout themselves. The example runner copies
`initSSLIP` and verifies after each example that all SSLIP functions, including
the example helpers, resolve inside its temporary checkout.

### Validation

The focused checks and full Ni/HCP examples run successfully from a clean copy
on MTEX 6.1.0 and 7.0.0. With rotation off and on, all eight sets of fitted
activities, residuals, solver flags, gradients, coordinates, slip tensors, and
optional rotations match their version's checkpoint-5 baseline exactly.
The 12 exported activity/rotation PNGs are pixel-identical.

A whole-object equality check differed for the saved MTEX 7 stress object in
one Ni rotation case. The physical stress matrix, rank, symmetry group,
coordinate-frame basis and names, and plotting convention match; all other
options match exactly. The numerical arrays and image pixels also match. The
comparison therefore checks these stress values explicitly rather than relying
on equality of the serialized MTEX object and its referenced objects.

An additional MTEX 7 run invoked each example by name from an unrelated working
folder, so MATLAB did not automatically switch into `examples/`. Both located
their inputs and helpers, retained their numerical results, and saved their MAT
files in that execution folder. A separate initializer check retained the
chosen MTEX and working directory and kept example utilities off the library
path until an example requested them.

All runs used disabled figure windows. The supplied Ni dataset remains tracked,
and local literature remains ignored. Source references and the adaptation
commit's `Co-authored-by` trailer retain Philipp's contribution credit.
Activity grouping and dominant-system ranking remain separate proposals.

## Reproducing the checks

See [tests/README.md](../tests/README.md). The example runner saves numeric
outputs independently of MTEX object serialization, including the solver exit
flags, along with an image of each example's final activity plot.
