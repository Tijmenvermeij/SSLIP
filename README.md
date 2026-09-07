## SSLIP

![SSLIP](https://ars.els-cdn.com/content/image/1-s2.0-S1359645422008795-ga1_lrg.jpg)

# Introduction to SSLIP

Slip System based Identification of Local Plasticity (**SSLIP**) is a methodology for automatic, point-by-point, identification of crystallographic slip system activity fields, performed by matching of Digital Image Correlation (DIC) displacement gradient fields to combinations of theoretical slip systems (resulting from e.g. EBSD). More details can be found in [**this paper**](https://doi.org/10.1016/j.actamat.2022.118502).

The **SSLIP** function library is written in [**MATLAB**](https://mathworks.com/products/matlab.html) and uses several functionalities of the MATLAB-based crystallographic toolbox [**MTEX**](https://mtex-toolbox.github.io). 

**Requires MTEX 6.1 or later. The regression checks and both examples have been tested with MTEX 6.1.0 and the official MTEX 7.0.0 release using MATLAB R2024b.** See the [validation results and numerical limitations](docs/mtex7-and-pr2-review.md).

Both MTEX versions use the same functions and options. Initialize your chosen
MTEX version in a fresh MATLAB session, then run the examples as usual. For a
quick compatibility check, add `tests/` to the MATLAB path and run
`runSSLIPChecks` ([background-run instructions](tests/README.md)).

It is important to use aligned EBSD/DIC data. See the following repository for an alignment framework: [**NanoMech_Alignment_Matlab**](https://github.com/Tijmenvermeij/NanoMech_Alignment_Matlab).

The **SSLIP** methodology and plotting functionalities are highlighted in a series of [example scripts](https://github.com/Tijmenvermeij/SSLIP#examples) that showcase how the functions work and what their output comprises.

Recently, this work was presented at the annual MTEX Workshop, see a recorded video [*here*](https://www.youtube.com/watch?v=xjNWsHeHnlA).

Please report any bugs you encounter. Instructions for the [regression checks and example runner](tests/README.md) are included in the repository.

# Authors
**SSLIP** has been conceptualized and created by [**Tijmen Vermeij**](https://www.tue.nl/en/research/researchers/tijmen-vermeij/), under supervision of **Johan Hoefnagels**, **Ron Peerlings** and **Marc Geers**.

**Philipp ([PhilKro](https://github.com/PhilKro))** contributed the single-slip
selection correction, explicit rotation support, and preprocessing separation adopted from
[PR #2](https://github.com/Tijmenvermeij/SSLIP/pull/2). The adapted commits retain
his co-author credit; the [integration notes](docs/mtex7-and-pr2-review.md)
identify his original commits and the corrections made during integration.

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

# How to cite SSLIP
If you have applied the SSLIP analyses to your research, please cite this open-access paper as your reference:
[**T. Vermeij, R.H.J. Peerlings, M.G.D. Geers, J.P.M. Hoefnagels, Automated identification of slip system activity fields from digital image correlation data, Acta Materialia: 243, 118502. (2022)**](https://doi.org/10.1016/j.actamat.2022.118502).
