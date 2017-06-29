## <a name="top">Meeting Notes

* [2017/06/29](#20170629)



-----------------------------------------------------------------
### Meeting on 2017/06/29

To cover at this and subsequent meetings:
- Goals of the project: near-, mid-, and long-term
- Structure of project
  - language(s): C++, python front end (c-extension? cython? need to choose how
    exactly we want to do this) [Sam]
  - third party libraries: deal-II for FEM; PyNE for data?
  - testing framework: google test (?), [Marissa]
  - integration framework: travis ci? batlab?  [Josh]
  - test problem library (benchmark suite)
  - build system: make
  - documentation: doxygen
  - style guide: TBD [Josh]
- Collaboration / interaction plan
- License [slaybaugh]
- Anything else?

#### Goals and ToDo Items

Goals: 
- research tool to investigate transport methods
- easy to get started as a new student
- easy to build
- easy to maintain
- easy to test
- easy to modify
- be modular; clear and simple API
- follows reproducible practices: git, tests, tags, documentation, style guide,
  build system, etc.
- Covariant reference frame

Include a tutorial

**near-term** (next 6 months)

Goal: structure in place; two methods available

Modules:
- Cartesian spatial solver (2D and/or 3D)
- access to nuclear data (via PyNE?)
- source specification
  - input 
  - internal representation
- geometry specification
- mesh
- boundary conditions (vacuum and incident)
- material assignment
- multiple energy groups


**mid-term** (out to 1.5 years)
- eigenvalue solver
- reflecting boundary conditions
- cylindrical


**long-term** (past 1.5 years)
- any additional boundary condition types
- multiple approaches for energy structure (?)
- more versatile geometry (?)
- spherical
- more methods of each category


[Index](#top)


