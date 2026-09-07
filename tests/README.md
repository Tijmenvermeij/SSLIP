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

Requires Optimization Toolbox and Parallel Computing Toolbox. The checks use
real `coneprog` solves but disable automatic pool creation temporarily, restoring
the previous setting afterwards. They cover absolute and relative residual
tolerances, the noise floor, signed and positive activity, skipped/invalid pixels,
infeasible fits, per-pixel solver exit flags, single-slip filtering and selection,
residual color scales, and all three trace plotting layouts.

## Full examples

The example runner also requires Image Processing Toolbox. It copies the code,
example scripts, and sample data into a temporary working folder inside the
specified output directory, then runs both complete examples. It saves numeric
results, including solver exit flags, and an image of each final activity plot.
The source checkout is not used for generated results.

```matlab
runSSLIPExamples('/path/to/results/mtex61');
```

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

To run the supplied scripts directly, use the `examples/` directory as the
working directory. They save result MAT files there.
