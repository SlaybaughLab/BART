[![Build Status](https://travis-ci.org/SlaybaughLab/BART.svg?branch=master)](https://travis-ci.org/SlaybaughLab/BART)

# Bay Area Radiation Transport (BART)

Software for investigating methods to solve the Boltzmann Transport equation.
Developed and maintained by 

- [Marissa Ramirez Zweiger](https://github.com/mzweig/)
- [Joshua Rehak](https://github.com/jsrehak/)
- [Rachel Slaybaugh](https://github.com/rachelslaybaugh)
- [Richard Vasques](https://github.com/ricvasques)
- [Weixiong Zheng](https://github.com/weixiong-zheng-berkeley/)

of the Nuclear Engineering Department at the University of California, Berkeley

## DG-EP

This is a project to do neutronics calculation using finite element methods.

The initial motivation is to explore the Discontinuous Galerkin (DG) method for the Even Parity
(EP) neutron transport equations in general dimensions (up to 3D). 
This is an initial setup and prototypical code with isotropic scattering.
The code is still in progress.

Up to now, the code is re-designed and gradually re-coded to take in different transport solvers in
an OOP way.

Currently, the project is named XTrans and will be discussed.

Currently, we are using [deal.II](http://www.dealii.org/) finite element
library. In the future, we are likely to convert to
[libMesh](http://libmesh.github.io/).

## Install and build
Please check [install_build.md](https://github.com/SlaybaughLab/BART/blob/master/install_build.md) for installation and building instructions.
