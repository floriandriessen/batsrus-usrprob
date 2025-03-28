# BATS-R-US user problems for Hot Jupiters

User problems to use with the [BATS-R-US](https://github.com/SWMFsoftware/BATSRUS)
code, version 9.20 (release 2018). Problems encompass single-fluid, multi-species radiation-(M)HD and multi-fluid radiation-(M)HD setups.

Additional code source files are included with modifications to output extra variables in the TecPlot I/O and to perform the plane-parallel ray-tracing calculation.

The problems only work with the 2018 release of BATS-R-US. Development had been started to port it to more recent version of the code, but this is unfinished work. Installation and compilation of the MHD source code can be read in the corresponding documentation. The latter documentation also explains how to setup a user problem as provided here.

## Useful webpages

- [2018 release](https://github.com/SWMFsoftware/BATSRUS/releases/tag/2018-06-24) binary of the code.

- [Nightly test page](http://herot.engin.umich.edu/~gtoth/) including links to various pieces of documentation.

## Problems

Each directory contains a setup to model the planetary atmosphere. A template input parameter file is also included. Further specifications on the included physics for each problem are given in the corresponding `ModUser.f90` file.

- Single-fluid, multi-species ideal MHD (`onefluid-multisp`). This code has been used for the study of [Presa, Driessen, & Vidotto (2024)](https://academic.oup.com/mnras/article/534/4/3622/7816390). It derives from the original code used in [Carolan et al. (2021)](https://academic.oup.com/mnras/article/508/4/6001/6395334), with the notable exception of fixing a bug of treating ion/neutral multi-species wrongly in the code. The bug was, however, minor since it only affected the flow near the planetary surface.

The remaining problems are all **unpublished work**. Most are still in a development phase:

- Two-fluid ideal MHD including electron pressure evolution (`twofluid-pe`). Evolves hydrogen ions and neutrals as separate fluids including the electron pressure as an additional hydro equation. It served to test whether collisional decoupling happened, especially of interest in the context of magnetic Hot Jupiters. Beware that this user problem is only tested during development.

- Four-fluid ideal MHD including electron pressure evolution (`fourfluid-pe-helium`). Evolves hydrogen and helium ions and neutrals as separate fluids including the electron pressure as an additional hydro equation. It served to test whether collisional decoupling and atmospheric fractionation happened. It was inspired by the work of [Xing, Yan, & Guo (2023)](https://iopscience.iop.org/article/10.3847/1538-4357/ace43f) that ran 1-D hydrogen-helium multi-fluid models with the PLUTO code. Beware that this user problem is only tested during development and some erratic behaviour developed. Particularly, the helium neutrals show inflow onto the planet from a large distance. Further investigation is warranted on what causes it.



