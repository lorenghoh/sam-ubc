# UBC Repository for System for Atmospheric Modelling (SAM)

This repository contains the source code for the UBC-maintained version 6.11.8 of the System for Atmospheric Modelling (SAM; Kharioudinov and Randall, 2003), which includes the (updated) `entrainment` module from Dawe and Austin (2011a) and a number of updates and improvements. The `entrainment` module itself also includes a few changes and fixes from the original paper, most notably a fix for calculations of advective fluxes, and saturation vapour pressure.

The main goal of this repository is to keep track of the changes made to the original model (Kharioudinov and Randall, 2003). Furthermore, future changes will be better documented and maintained. This is useful as we plan to add a number of helper scripts that can streamline the model use, which will be frequently updated.

This version also includes a newly constructed build system based on CMake (refer to `README_CMAKE.md` for detailed build instructions ). 


## References

Dawe, J. T., & Austin, P. H. (2011). Interpolation of LES cloud surfaces for use in direct calculations of entrainment and detrainment. *Monthly weather review*, **139(2)**, 444-456, **https://doi.org/10.1175/2010MWR3473.1**

Khairoutdinov, M. F., & Randall, D. A. (2003). Cloud Resolving Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, Uncertainties, and Sensitivities. *Journal of the Atmospheric Sciences*, **60(4)**, 607–625. **https://doi.org/10.1175/1520-0469(2003)060%3C);0607:crmota>2.0.co;2**
