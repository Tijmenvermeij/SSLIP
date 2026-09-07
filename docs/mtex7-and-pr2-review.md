# MTEX 7 validation and Philipp's PR #2

Date: 7 September 2026. Cleanup checkpoint: `f0ef743`.

## Decision

Keep the existing public interface and code organization for this compatibility
update. The regression checks and both examples run under MTEX 6.1.0 and the
[official MTEX 7.0.0 release](https://github.com/mtex-toolbox/mtex/releases/tag/mtex-7.0.0).
No broad MTEX API rewrite was needed for these paths.

Port Philipp's `singleSlipPerPixel` correction. Defer merging the full PR until
the interface, indexing, and rotation issues below are resolved. The review
branch `codex/philipp-pr2-review` preserves his original revision.

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

## PR #2 review

Reviewed [PR #2](https://github.com/Tijmenvermeij/SSLIP/pull/2),
"Implemented rotation, straightened out architecture", at
`222a9de2448eb7e152a422d1b3bf0d335fd2bb06`.
Its merge base with current main is `1e2faa0` from 3 December 2025, so integration
must preserve subsequent main changes, including the newer example dataset.

### Confirmed issues to fix before merging

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
- Explicit rotation control is preferable to overloading `CRSS == 0`.
  The PR penalizes rotation coefficients in the same L1 objective as slip, and
  still skips low-effective-strain pixels, including pure rotation when the
  threshold is positive. Those choices need to be specified and tested before
  calling this a general rotation correction.
- Preserve the current sample dataset and citation information during a future
  merge. The PR removes the sample MAT file and removes the paper citation
  header from the combined solver.

## Reproducing the checks

See [tests/README.md](../tests/README.md). The example runner saves numeric
outputs independently of MTEX object serialization, including the solver exit
flags, along with an image of each example's final activity plot.
