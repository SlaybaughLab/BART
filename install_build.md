# Installations
## Packages requirement on BART
The C++ core is based upon the third-party library [deal.II](http://www.dealii.org/). Besides, PETSc is required to provided MPI interface and corresponding objects such as sparsematrix etc.

Before installation process, it is required to have CMake. If installation is on Mac OSX, xcode and commend line tool is also required. Checkout deal.II [wiki](https://github.com/dealii/dealii/wiki/MacOSX) for instructions.

## Mac OSX through .dmg
deal.II installation package can be downloaded [here](https://www.dealii.org/download.html) on the upper right corner. Checkout the Section "Detailed instructions" of deal.II [wiki page](https://github.com/dealii/dealii/wiki/MacOSX) for installation steps.

This installation contains every aspect of deal.II including PETSc, Trilinos, etc.
## Linux/Unix and Mac OSX through Homebrew/Linuxbrew
Installation on Linux/Unix and OSX can be done through Homebrew Linuxbrew as well. Checkout [here](https://github.com/dealii/dealii/wiki/deal.II-on-Homebrew---Linuxbrew) for instructions.

## Installation through candi
Installation can also be done through following instructions in [candi](https://github.com/koecher/candi).

## Build from scratch
This is to build deal.II with PETSc from scratch by oneself, which gives the most flexibility but is not recommended. If interested, please check [here](https://www.dealii.org/developer/readme.html).
# Build BART
Up to 2017-07-13, one could only build the core of BART, which is [XTrans](https://github.com/weixiong-zheng-berkeley/XTrans). In the future this will be updated.

To build the core, one needs to add two lines to the .bash type file to tell where your deal.II build is:

`export DEAL_II_CONF_SILENT=ON`

`. your_path_to_dealii`

And then source it. In the source code directory (where there is the CMakeLists.txt), type:

`cmake .`

`make`

Notice that in the future for running research problems, rather than development, replace `make` with:

`make release`

so that the code will run in release mode, which can be several times faster.
