# SAM_ENT

A Git repository for the System for Atmospheric Modelling ([SAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html)) with the entrainment module developed by Dawe and Austin (2011). It is a large-eddy simulation (LES) model developed by Khairoutdinov and Randall (2003), used to perform realistic simulations of the cloudy atmosphere.

We have applied a few modifications to SAM. We added the *entrainment* module for the direct calculation of entrainment and detrainment rates using the tetrahedral interpolation scheme (Dawe and Austin, 2011), which is used for a precise quantification of entrainment and detrainment rates during moist convection. The entrainment module itself also includes a few changes and fixes from the original paper, most notably a fix for calculations of advective fluxes, and saturation vapour pressure.

We have also added CmakeList.txt file for a faster and more stable build process. Take a look at README_CMAKE.md doc for a detailed description.


Dawe, J. T., & Austin, P. H. (2011). Interpolation of LES cloud surfaces for use in direct calculations of entrainment and detrainment. *Monthly weather review*, **139(2)**, 444-456.

Khairoutdinov, M. F., & Randall, D. A. (2003). Cloud resolving modeling of the ARM summer 1997 IOP: Model formulation, results, uncertainties, and sensitivities. *Journal of the Atmospheric Sciences*, **60(4)**, 607-625.
