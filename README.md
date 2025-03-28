# BATS-R-US user problems for Hot Jupiters

User problems to use with the [BATS-R-US](https://github.com/SWMFsoftware/BATSRUS)
code. Problems encompass single-fluid, multi-species radiation-(M)HD and multi-fluid radiation-(M)HD setups.

The problems work with the 2024 release of BATS-R-US (at the time of writing any recent version works). Installation and compilation of the MHD source code can be read in the corresponding documentation. The latter documentation also explains how to setup a user problem as provided here.

## Useful webpages

- [Nightly test page](http://herot.engin.umich.edu/~gtoth/) including links to various pieces of documentation.

## Problems

Each directory contains a setup to model the planetary atmosphere. A template input parameter file is also included. Further specifications on the included physics for each problem are given in the corresponding `ModUser.f90` file.

- Single-fluid, multi-species ideal MHD (`onefluid-multisp`). This code has been used for the study of [Presa, Driessen, & Vidotto (2024)](https://academic.oup.com/mnras/article/534/4/3622/7816390). It derives from the original code used in [Carolan et al. (2021)](https://academic.oup.com/mnras/article/508/4/6001/6395334), with the notable exception of fixing a bug of treating ion/neutral multi-species wrongly in the code. The bug was, however, minor since it only affected the flow near the planetary surface.



