| `master` | `dev` |
|----------|-------|
|[![Build Status](https://travis-ci.org/SlaybaughLab/BART.svg?branch=master)](https://travis-ci.org/SlaybaughLab/BART)[![Codecov](https://codecov.io/gh/SlaybaughLab/BART/branch/master/graph/badge.svg)](https://codecov.io/gh/SlaybaughLab/BART/branch/master)|[![Build Status](https://travis-ci.org/SlaybaughLab/BART.svg?branch=dev)](https://travis-ci.org/SlaybaughLab/BART)[![Codecov](https://codecov.io/gh/SlaybaughLab/BART/branch/dev/graph/badge.svg)](https://codecov.io/gh/SlaybaughLab/BART/branch/dev)|

# Bay Area Radiation Transport (BART)

## What is BART for?

Bay Area Radiation Transport (BART), is a C++-based research-purpose
code for parallel radiation transport computation in nuclear reactor
applications. BART is under active developed by the Computational
Neutronics group in Nuclear Engineering at University of California,
Berkeley.

### Documentation

BART documentation is generated using [doxygen](http://www.stack.nl/~dimitri/doxygen/).

### Test-driven development
BART uses the following applications and libraries for testing:
- [Travis CI](https://travis-ci.org) for continuous integration
- [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest) for unit testings requiring MPI
- [Google Test](https://github.com/google/googletest) for all other
  serial unit testings.
- [Codecov](https://codecov.io/) for code coverage of unit tests.

## Rationale
### Finite element method based code
BART is based off the general purpose finite element method [deal.II](http://www.dealii.org/). It aims to solve first and second-order forms of linear Boltzmann equation for nuclear reactor applications using continuous/discontinuous finite element methods for spatial discretization with existing/developing acceleration methods. BART uses discrete ordinates method for angular discretization. 

### Parallelism in meshing and linear algebra
BART uses MPI for parallelism and is designed for computation on distributed memory system:
- By utilizing distributed triangulation enabled by [p4est](https://www.mcs.anl.gov/petsc/) library wrapped in deal.II, BART can automatically partition the mesh by however many number of processors requested and distribute the triangulation onto different processors.
- BART heavily depends on [PETSc](https://www.mcs.anl.gov/petsc/) by utilizing deal.II wrappers of PETSc data structure. Therefore, all the parallel-supported functionalities in PETSc (if wrapped by deal.II) can be invoked by BART. This includes parallel sparse matrix, parallel vectors and parallel preconditioners/algebraic solvers.

### Formulations

BART supports the Diffusion Equation in 1/2/3D, and second-order forms of the transport equation (such as the Self-Adjoint Angular Flux and Even Parity) will be implemented.

### Acelleration Methods

One of the major design goals of BART is to provide a framework for testing acceleration methods. There are no methods implemented in the current version, but multiple methods are planned to be implemented, including:

- Nonlinear diffusion acceleration (NDA)
- Two-grid acceleration (TG)

### Benchmarks

Benchmarks from Sood (1999) are provided in the `benchmarks` folder for validation of the code.

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

## Developers
The development work is led by [Dr. Rachel Slaybaugh](https://github.com/rachelslaybaugh). Graduate students actively involved in the development include:
- [Joshua Rehak](https://github.com/jsrehak/)

Previous developers include 
- [Dr. Weixiong Zheng](https://github.com/weixiong-zheng-berkeley/).
- [Marissa Ramirez Zweiger](https://github.com/mzweig/)
- [Alexander Blank](https://github.com/AlexanderBlank)
