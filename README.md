# SAM_ENT

A Git repository for a modified version of the System for Atmospheric Modelling ([SAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html)). It is a large-eddy simulation (LES) model developed by Khairoutdinov and Randall (2003), used to perform realistic simulations of the cloudy atmosphere.

We have applied a few modifications to SAM. We added the direct entrainment calculation using tetrahedral interpolation (Dawe and Austin, 2010), which is used for a precise quantification of entrainment and detrainment rates during moist convection. The entrainment module itself also includes a few changes and fixes from the original paper, most notably a fix for calculations of advective fluxes, and saturation vapour pressure.

We have also added CmakeList.txt file for a faster and more stable build process. Take a look at README_CMAKE.md doc for a detailed description.
