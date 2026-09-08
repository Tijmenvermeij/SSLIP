# Regression and compatibility checks

Use a fresh MATLAB session with the MTEX version to be tested initialized.
Avoid switching MTEX versions in an already-running session. The plotting
checks and supplied examples close figures.

For background runs, start MATLAB with `-noFigureWindows -nosplash` so no plot
windows appear, even briefly. MTEX 6.1 explicitly makes figures visible, so
`DefaultFigureVisible = 'off'` alone is insufficient. Plot files can still be
exported with figure windows disabled. With the desired MTEX initialized by
your MATLAB startup script:

```sh
matlab -noFigureWindows -nosplash -batch "addpath('/path/to/SSLIP/tests'); runSSLIPChecks;"
```

```matlab
addpath('/path/to/mtex');
startup_mtex('noMenu');
addpath('/path/to/SSLIP/tests');
runSSLIPChecks;
```

The checks add `src` and its subfolders with `addpath(genpath(...))`
automatically. No separate initialization function is required.

Requires Optimization, Parallel Computing, and Image Processing Toolboxes. The checks use
real `coneprog` solves but disable automatic pool creation temporarily, restoring
the previous setting afterwards. They cover absolute and relative residual
tolerances, the noise floor, signed and positive activity, skipped/invalid pixels,
infeasible fits, per-pixel solver exit flags, single-slip filtering and selection,
residual color scales, and all three trace plotting layouts.

Stress-alignment checks cover a general stress tensor, equivalent uniaxial
vector/tensor inputs, reordered subsets, the full returned fitting basis,
unchanged CRSS, normalized activities, rotation, and reconstruction for methods
1 and 2. They also check zero resolved shear, zero stress, unchanged signed
fits and method 3, invalid inputs, and the warning when stress is omitted.

Preprocessing checks cover affine gradients on a shifted rectangular grid, shuffled
EBSD input, coarse-graining, zero and missing displacement values with filtering
off and on, and rejection of mismatched displacement dimensions.

Direct-gradient checks cover structured input, shuffled point order, preserved
zero and missing components, incomplete/mismatched/ambiguous data, unsupported
preprocessing, rotation fitting, removal of stale result fields, and five-panel
gradient-only plots including all-zero data. Structured displacement input is
also checked against the original five-input call.

Single-slip checks cover reordered subsets, the reported positive-constraint
override, tied residuals, normalized coefficients, matrix-shaped gradient input,
and strict rejection at the residual threshold. They retain method 3's existing
behavior for low-strain and missing pixels; method 1's `minEeff` cutoff is not
applied to the single-slip fit.

Separate plotting checks cover custom deformation limits, logarithmic shared
residual scales, use of the supplied residual when stored data differs, rotation
display in degrees without changing the saved radians, missing rotation data,
and the existing PNG and additional JPEG export names. These use temporary
output folders and restore the working directory afterwards.

The same checks also exercise Philipp's optional rotation extension: positive
and negative rotation with either slip-sign constraint, radians after in-plane
normalization, pure slip, pure rotation, mixtures, a noisy mixture, skipped and
invalid pixels, selected-system plotting, and rejection by unsupported methods.
These synthetic fits use independent slip and rotation bases so known
amplitudes can be recovered; they do not imply uniqueness for arbitrary sets
of crystallographic slip systems.

## Full examples

The example runner also requires Image Processing Toolbox. It copies the code,
example scripts and helpers, and sample data into a temporary working folder
inside the specified output directory, then runs both complete examples. It saves numeric
results, including solver exit flags, and an image of each final activity plot.
After each example, it verifies that every SSLIP function and example helper
resolves inside the temporary copy. This catches missing paths or accidental
dependence on another checkout. The source checkout is not used for generated
results.

```matlab
runSSLIPExamples('/path/to/results/mtex61');
```

The examples explicitly set `IDoptions.enableRotation = 0`. To test their
rotation option and rotation plots, pass `true`; only the temporary copies are
changed. Use a separate output directory to retain the default results:

```matlab
runSSLIPExamples('/path/to/results/mtex61-rotation', true);
```

This also saves the fitted rotation in the numeric MAT files (radians) and
exports each example's rotation plot (degrees). The physical slip plots keep
their original system order and zero-activity handling.

Repeat in a fresh session with MTEX 7 initialized and a different output directory:

```matlab
runSSLIPChecks;
runSSLIPExamples('/path/to/results/mtex70');
comparison = compareSSLIPExamples('/path/to/results/mtex61', ...
                                 '/path/to/results/mtex70');
```

The comparison checks grid and slip-system alignment and reports solver-status,
activity, residual, reconstructed-gradient, and objective differences. It does
not assume identical slip amplitudes or silently accept differing solver status.
See the [validation report](../docs/mtex7-and-pr2-review.md) for the measured
round-off sensitivity of the example problems.

To run a supplied script directly, use its full path, for example:

```matlab
run('/path/to/SSLIP/examples/NiSuperAlloyExperiment.m');
```

The script locates its own code and data. MATLAB's `run` executes it in
`examples/`, so its result MAT file is saved there.
