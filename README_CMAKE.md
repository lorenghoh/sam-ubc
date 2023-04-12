Getting started with SAM_CMAKE.

The following documentation provides best practices for running SAM
on a high-performance cluster. The following steps can be automated
using the scripts in cmake_snippets repository:

> https://github.com/lorenghoh/sam_snippets

but you will still be responsible for creating a case folder containing
all the necessary files/folders/symlinks. 


 1. Set up an out-of-source compile location
    * For Optimum, the case directory must be off /scratch and /data
    e.g. mkdir ~/scratch/BOMEX/

 2. Edit SRC/domain.f90 to set domain size and number of sub-domains

 3. Edit CMakeLists.txt and set the corresponding netcdf and openmpi
    library links (Intel compilers recommended)
 
 4. Use cmake_snippets from the repository
        i.e. https://github.com/lorenghoh/sam_snippets

    Or, compile directly from source
        e.g. cmake /home/loh/SAM

 5. The cmake script should have automatically linked all the necessary
files, but you are responsible for copying and modyfying case-related
files (along with SAM/CaseName)

 6. Create output directories:
    mkdir -p OUT_STAT
    mkdir -p OUT_3D
    mkdir -p RESTART

    For Optimum, OUT_3D should most likely be in /scratch, in which 
    case one would instead create a symlink to the scratch space:
        e.g. ln -s /scratch/BOMEX OUT_3D

    and make sure the Lustre striping is turned off.
        e.g. lfs setstripe -c 1 /scratch/BOMEX

 7. Run the model. For example, given a 16-task job,
    mpiexec -n 16 SAM_CMAKE
