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

## How do we manage the development?

### Documentation

BART documentation is generated using [doxygen](http://www.stack.nl/~dimitri/doxygen/).

### Test-driven development
BART uses the following applications and libraries for testing:
- [Travis CI](https://travis-ci.org) for continuous integration
- [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest) for unit testings requiring MPI
- [Google Test](https://github.com/google/googletest) for all other
  serial unit testings.
- [Codecov](https://codecov.io/) for code coverage of serial tests.

<!-- ### Agile management -->
<!-- We are gradually immersing ourselves in the principle of agile management using [Jira](https://www.atlassian.com/software/jira?aceid=&adposition=1t1&adgroup=9124375582&campaign=189421462&creative=256725234926&device=c&keyword=jira&matchtype=e&network=g&placement=&ds_kids=p19481846873&gclid=Cj0KCQjwkpfWBRDZARIsAAfeXarkxD2j0JPwTTaH07dxEy8nVbZgK7U_Uj8hDx7j2uyUBXl29zrtoQQaAshhEALw_wcB&gclsrc=aw.ds) to improve our BART development tracking. -->

## Developers
The development work is led by [Dr. Rachel Slaybaugh](https://github.com/rachelslaybaugh). Graduate students actively involved in the development include:
- [Joshua Rehak](https://github.com/jsrehak/)

Previous developers include [Dr. Weixiong Zheng](https://github.com/weixiong-zheng-berkeley/).

## What is the rationale behind BART?
### BART is a finite element method based code
BART is based off the general purpose finite elment [deal.II](http://www.dealii.org/). It aims to solve first and second-order forms of linear Boltzmann equation for nuclear reactor applications using continuous/discontinuous finite element methods for spatial discretization with existing/developing acceleration methods. BART uses discrete ordinates method for angular discretization. 

### Parallelism in meshing and linear algebra
BART uses MPI for parallelism and is designed for computation on distributed memory system:
- By utilizing distributed triangulation enabled by [p4est](https://www.mcs.anl.gov/petsc/) library wrapped in deal.II, BART can automatically partition the mesh by however many number of processors requested and distribute the triangulation onto different processors.
- BART heavily depends on [PETSc](https://www.mcs.anl.gov/petsc/) by utilizing deal.II wrappers of PETSc data structure. Therefore, all the parallel-supported functionalities in PETSc (if wrapped by deal.II) can be invoked by BART. This includes parallel sparse matrix, parallel vectors and parallel preconditioners/algebraic solvers.

### Supported Transport Equation forms

Bart supports the following forms of the transport equation:

- Self Adjoint Angular Flux
- Even Parity

More forms of the transport equation and accelleration methods are an
area of active development.

### Meshing capability
Supported meshes in BART include:
- Hyper-rectangular mesh in 1/2/3D;
- Fuel pin-resolved curvilinear mesh in 2D;
- Fuel pin-resolved curvilinear mesh in 3D based on extrusion.

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
