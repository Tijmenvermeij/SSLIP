## SSLIP

![SSLIP](https://ars.els-cdn.com/content/image/1-s2.0-S1359645422008795-ga1_lrg.jpg)

# Introduction to SSLIP

Slip System based Identification of Local Plasticity (**SSLIP**) is a methodology for automatic, point-by-point, identification of crystallographic slip system activity fields, performed by matching of Digital Image Correlation (DIC) displacement gradient fields to combinations of theoretical slip systems (resulting from e.g. EBSD). More details can be found in [**this paper**](https://doi.org/10.1016/j.actamat.2022.118502).

The **SSLIP** function library is written in [**MATLAB**](https://mathworks.com/products/matlab.html) and uses several functionalities of the MATLAB-based crystallographic toolbox [**MTEX**](https://mtex-toolbox.github.io). 

**Requires MTEX 6.1 or later. The regression checks and both examples have been tested with MTEX 6.1.0 and the official MTEX 7.0.0 release using MATLAB R2024b.** See the [validation results and numerical limitations](docs/mtex7-and-pr2-review.md).

Both MTEX versions use the same functions and options. Initialize your chosen
MTEX version in a fresh MATLAB session, then initialize SSLIP as shown below.
For a quick compatibility check, add `tests/` to the MATLAB path and run
`runSSLIPChecks` ([background-run instructions](tests/README.md)).

It is important to use aligned EBSD/DIC data. See the following repository for an alignment framework: [**NanoMech_Alignment_Matlab**](https://github.com/Tijmenvermeij/NanoMech_Alignment_Matlab).

The **SSLIP** methodology and plotting functionalities are highlighted in a series of [example scripts](https://github.com/Tijmenvermeij/SSLIP#examples) that showcase how the functions work and what their output comprises.

Recently, this work was presented at the annual MTEX Workshop, see a recorded video [*here*](https://www.youtube.com/watch?v=xjNWsHeHnlA).

Please report any bugs you encounter. Instructions for the [regression checks and example runner](tests/README.md) are included in the repository.

# Authors
**SSLIP** has been conceptualized and created by [**Tijmen Vermeij**](https://www.tue.nl/en/research/researchers/tijmen-vermeij/), under supervision of **Johan Hoefnagels**, **Ron Peerlings** and **Marc Geers**.

**Philipp ([PhilKro](https://github.com/PhilKro))** contributed the single-slip
selection correction, explicit rotation support, and separation of preprocessing,
single-slip solving, and the small plotting routines adopted from
[PR #2](https://github.com/Tijmenvermeij/SSLIP/pull/2). The adapted commits retain
his co-author credit; the [integration notes](docs/mtex7-and-pr2-review.md)
identify his original commits and the corrections made during integration.
The structured displacement/gradient input, automatic stress alignment, folder
layout, and `initSSLIP` initializer are also adapted from his PR.

# Getting started

Start your chosen MTEX version, then add the SSLIP repository root and run
Philipp's initializer:

```matlab
addpath('/path/to/SSLIP');
initSSLIP;
```

This adds `src`, `src/plotting`, and `src/utils` for the current session.
Replace old `addpath('/path/to/SSLIP/src')` setup with the two lines above;
adding only `src` no longer includes the moved plotting and utility functions.
The initializer leaves MTEX selection, the working directory, and the saved
MATLAB path to the caller.

The adopted folder layout is:

```text
initSSLIP.m       Library initialization
src/             SSLIP entry point, preprocessing, and solvers
src/plotting/    Activity, deformation, residual, and rotation plots
src/utils/       Shared numerical and grid helpers
examples/        Ni and virtual HCP example scripts
examples/utils/  Synthetic slip-step generation helpers
data/            Supplied Ni dataset
tests/           Regression checks and example runner
```

# Examples

After starting MTEX, run either example by its full path:

```matlab
run('/path/to/SSLIP/examples/NiSuperAlloyExperiment.m');
run('/path/to/SSLIP/examples/virtualExperimentHCP.m');
```

The scripts initialize SSLIP and add their own example helpers. They locate
code and data relative to the script file, so input lookup does not depend
on the current working directory. Both scripts clear workspace variables and
close figures before running.

Results are saved in the directory where the script executes. MATLAB's `run`
command above executes in `examples/` and saves the result MAT files there.
For isolated copies and exported comparison figures, use the
[example runner](tests/README.md#full-examples).

# Optional rotation correction

Set `opt.enableRotation = 1` with `opt.IDMethod = 1` to use Philipp's rotation
extension in the existing `SSLIP(ebsd,U,V,sSLocal,opt)` call. Rotation is returned
separately in `ebsdID.prop.rotationIDcor`, in radians, while slip activities keep
their existing shape and system order. Rotation remains signed even with
`opt.posConstr = 1`. Methods 2 and 3 reject this option.

Both supplied examples expose this as `IDoptions.enableRotation` (default 0).
Set it to 1 to add a rotation plot in degrees; the saved numerical angle remains
in radians. Their rotation-enabled output names include `_rotation` to keep
them separate from the default results. In the virtual HCP example, batch 2
selects the actual generated slip systems `[3 7 19]`.

The existing `minEeff` cutoff still applies; use `opt.minEeff = 0` when pure
rotation should also be fitted. This uses the small-angle approximation from
[Vermeij et al., Strain (2025), Eq. (11)](https://doi.org/10.1111/str.70000).
It retains the PR's L1 penalty on rotation and does not implement the full
Radon and slip-pair selection workflow from that paper. See the
[formulation and normalization details](docs/mtex7-and-pr2-review.md#selected-rotation-integration).

# Slip signs for positive activity

For methods 1 and 2, set `opt.posConstr = 1` and supply `opt.stress` to use
Philipp's automatic alignment. SSLIP flips each Burgers vector when its
resolved shear stress is negative. Supply one `stressTensor` in specimen
coordinates, or one `vector3d` for a uniaxial tension direction:

```matlab
opt.posConstr = 1;
opt.stress = stressTensor.uniaxial(xvector);
[ebsdID,optOut,sSLocal] = SSLIP(ebsd,U,V,sSLocal,opt);
save(optOut.plotname,'ebsdID','optOut','sSLocal');
```

The third output is the **full list of systems after alignment**. Save and use
it for reconstruction and plotting; activity row `k` belongs to
`sSLocal(optOut.NoSs(k))`. This also works with the structured input below.
Both examples now capture and save this output and specify their tensile load.

Only the sign of the resolved shear is used; stress magnitude does not weight
the fit. Systems with exactly zero resolved shear retain their direction.
Without `opt.stress`, a positive-constrained fit warns and uses the supplied
systems. Signed fits (`posConstr = 0`) retain the supplied directions; method 3
also retains them because it always switches the positive constraint off.
Rotation remains signed. Use a negative stress tensor to represent compression;
negating a tension direction still represents the same uniaxial tension.

# Input from displacement gradients

The original `SSLIP(ebsd,U,V,sSLocal,opt)` call is retained. Philipp's structured
input form is also supported, with flat options:

```matlab
deformationData = struct('Hxx',Hxx,'Hxy',Hxy,'Hyx',Hyx,'Hyy',Hyy);
opt.filterSize = 0;
opt.coarsegrain = 1;
[ebsdID,optOut,sSLocal] = SSLIP(ebsd,deformationData,sSLocal,opt);
```

Supply displacement gradients: `Hxx = dU/dx`, `Hxy = dU/dy`, `Hyx = dV/dx`,
and `Hyy = dV/dy`. Each array must match the EBSD grid's size and point order.
Values must be real floating-point numbers; zeros are retained and NaN denotes
missing data. Single-precision inputs are converted to double for fitting.

These gradients must be ready for fitting. The gradient-input defaults are
`filterSize = 0` and `coarsegrain = 1`; requesting further filtering or
coarse-graining raises an error. Supply the corresponding processed EBSD grid
if your gradients have already been coarse-grained. This avoids guessing the
grid or applying preprocessing twice.

Alternatively, `deformationData = struct('U',U,'V',V)` uses the existing
displacement preprocessing and defaults. Mixing displacements and gradients
in the same input is rejected. With gradient-only input, the output contains
no `U` or `V`, and the deformation figure shows the five available fields.

# Replot saved results

The `ebsdID` and `optOut` variables saved by the examples can now be used to
redraw individual figures without rerunning identification:

```matlab
plotSSLIP_DeformationFields(ebsdID,optOut);
plotSSLIP_Residual(ebsdID,optOut);
rotationFigure = plotSSLIP_Rotation(ebsdID,optOut);
```

These retain the existing figure scales and export filenames. Rotation is
displayed in degrees; the stored values stay in radians. If no rotation field
is present, the rotation plotter returns an empty handle and creates no figure.
The existing `plotSSLIP` call still draws the activity maps and optional residual.

# How to cite SSLIP
If you have applied the SSLIP analyses to your research, please cite this open-access paper as your reference:
[**T. Vermeij, R.H.J. Peerlings, M.G.D. Geers, J.P.M. Hoefnagels, Automated identification of slip system activity fields from digital image correlation data, Acta Materialia: 243, 118502. (2022)**](https://doi.org/10.1016/j.actamat.2022.118502).
