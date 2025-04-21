## SSLIP

![SSLIP](https://ars.els-cdn.com/content/image/1-s2.0-S1359645422008795-ga1_lrg.jpg)

# Introduction to SSLIP

Slip System based Identification of Local Plasticity (**SSLIP**) is a methodology for automatic, point-by-point, identification of crystallographic slip system activity fields, performed by matching of Digital Image Correlation (DIC) displacement gradient fields to combinations of theoretical slip systems (resulting from e.g. EBSD). More details can be found in [**this paper**](https://doi.org/10.1016/j.actamat.2022.118502).

The **SSLIP** function library is written in [**MATLAB**](https://mathworks.com/products/matlab.html) and uses several functionalities of the MATLAB-based crystallographic toolbox [**MTEX**](https://mtex-toolbox.github.io). 

**WARNING: Make sure to use Mtex 5.11 or lower, Mtex 6.0 and above does not work with SSLIP**

It is important to use aligned EBSD/DIC data. See the following repository for an alignment framework: [**NanoMech_Alignment_Matlab**](https://github.com/Tijmenvermeij/NanoMech_Alignment_Matlab).

The **SSLIP** methodology and plotting functionalities are highlighted in a series of [example scripts](https://github.com/Tijmenvermeij/SSLIP#examples) that showcase how the functions work and what their output comprises.

Recently, this work was presented at the annual MTEX Workshop, see a recorded video [*here*](https://www.youtube.com/watch?v=xjNWsHeHnlA).

Please report any bugs you encounter.

# Authors
**SSLIP** has been conceptualized and created by [**Tijmen Vermeij**](https://www.tue.nl/en/research/researchers/tijmen-vermeij/), under supervision of **Johan Hoefnagels**, **Ron Peerlings** and **Marc Geers**.

# How to cite SSLIP
If you have applied the SSLIP analyses to your research, please cite this open-access paper as your reference:
[**T. Vermeij, R.H.J. Peerlings, M.G.D. Geers, J.P.M. Hoefnagels, Automated identification of slip system activity fields from digital image correlation data, Acta Materialia: 243, 118502. (2022)**](https://doi.org/10.1016/j.actamat.2022.118502).

