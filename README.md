[![Build Status](https://travis-ci.org/SlaybaughLab/BART.svg?branch=master)](https://travis-ci.org/SlaybaughLab/BART)

# Bay Area Radiation Transport (BART)
## What is BART for?
Bay Area Radiation Transport, AKA BART, is a C++ based research-purpose code for parallel radiation transport computation in nuclear reactor applications. BART is being actively developed at Computational Neutronics group in Nuclear Engineering at University of California, Berkeley.

## How do we manage the development?
### Documentation
BART documents everything implemented for functionality using [doxygen](http://www.stack.nl/~dimitri/doxygen/). We believe that we and future developers will be thankful to ourselves for documenting today at some time in the future. 

### Test-driven development
BART in the [restart run](https://github.com/SlaybaughLab/BART/tree/restart) is driven by the principle of testing everything necessary. Specifically, we use:
- [Travis CI](https://circleci.com/?utm_source=gnb&utm_medium=SEM&utm_campaign=SEM-gnb-400-Eng-ni&utm_content=SEM-gnb-400-Eng-ni-Travis_CI&gclid=Cj0KCQjwkpfWBRDZARIsAAfeXaqgCfU-RRPCdyxlvRTTGYw2qT31LlLfwIw_1OXVUw5gYvJSZhcHG9saAm_nEALw_wcB) for continuous integration;
- [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest) for unit testings requiring MPI;
- [Google Test](https://github.com/google/googletest) for all other serial unit testings.

## Who do the hard work?
The development work is led by -[Dr. Weixiong Zheng](https://github.com/weixiong-zheng-berkeley/) and -[Dr. Rachel Slaybaugh](https://github.com/rachelslaybaugh). Graduate students actively involved in the development include:
- [Joshua Rehak](https://github.com/jsrehak/)
- [Marissa Ramirez Zweiger](https://github.com/mzweig/)
In the long term, BART will welcome anyone who has a good idea to participating.

## What is the rationale behind BART?
### BART is a finite element method based code
BART is based off the general purpose finite elment [deal.II](http://www.dealii.org/). It aims to solve first and second-order forms of linear Boltzmann equation for nuclear reactor applications using continuous/discontinuous finite element methods for spatial discretization with existing/developing acceleration methods. BART uses discrete ordinates method for angular discretization. 

### Parallelism in meshing and linear algebra
BART is using MPI for parallelism. BART is designed for computation on distributed memory system:
- By utilizing distributed triangulation enabled by [p4est](https://www.mcs.anl.gov/petsc/) library wrapped in deal.II, BART can automatically partition the mesh by however many number of processors requested and distribute the triangulation onto different processors.
- BART heavily depends on [PETSc](https://www.mcs.anl.gov/petsc/) by utilizing deal.II wrappers of PETSc data structure. Therefore, all the parallel-supported functionalities in PETSc (if wrapped by deal.II) can be invoked by BART. This includes parallel sparse matrix, parallel vectors and parallel preconditioners/algebraic solvers.

### General dimensionality
BART was initially implemented for 2D and 3D for parallel computation. In the -[restart branch](https://github.com/SlaybaughLab/BART/tree/restart), the functionality is generalized to 1D with serial settings.

### Meshing capability
Originally, BART was implemented for homogenized mesh using rectangle mesh in 2D and regular cubois mesh in 3D. In the -[restart branch](https://github.com/SlaybaughLab/BART/tree/restart), some thrilling new features are implemented. Overall, we have:
- Hyper-rectangular mesh in 1/2/3D;
- Fuel pin-resolved curvilinear mesh in 2D;
- Fuel pin-resolved curvilinear mesh in 3D based on extrusion.

The thrilling part is that pin-resolved meshing does not require third-party library, e.g. Cubit and GMSH but rather similar to regular homogenized mesh.

Part of the work also contributes to development version of [deal.II](http://www.dealii.org/).

# Install and build
Please check [install_build.md](https://github.com/SlaybaughLab/BART/blob/master/install_build.md) for installation and building instructions.

# More to read
- [deal.II](http://www.dealii.org/): general purpose finite element library
- [PETSc](https://www.mcs.anl.gov/petsc/): high-performance parallel linear algebra library
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/): documentation.
- [p4est](http://www.p4est.org/): distributed meshing in multi-D.
- [Google Style Guide](https://google.github.io/styleguide/cppguide.html): consistent code convention.
- [Google Test](https://github.com/google/googletest): efficient unit testing tools.
- [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest): unit testing tools  with Google Test.
