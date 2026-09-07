# Regression checks

Initialize the MTEX version to be tested, add this directory to the MATLAB
path, and run `runSSLIPChecks`. Use a separate MATLAB session because the
plotting checks close figures.

```matlab
addpath('/path/to/mtex');
startup_mtex('noMenu');
addpath('/path/to/SSLIP/tests');
runSSLIPChecks;
```

Requires Optimization Toolbox and Parallel Computing Toolbox. The checks use
real `coneprog` solves but disable automatic pool creation for this session
only. They cover absolute and relative residual tolerances, the noise floor,
positive and signed slip activity, invalid and low-strain pixels, single-slip
filtering, residual color scales, and all three trace plotting layouts.

For a compatibility check, also run both scripts in `examples/` from that
directory. The example scripts save result MAT files to their working directory.
