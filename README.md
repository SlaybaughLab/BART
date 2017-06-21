# DG-EP

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
