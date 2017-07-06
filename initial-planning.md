## <a name="top">Meeting Notes

* [Goals](#goals)
* [Structure](#structure)
* [To Do](#todo)


-----------------------------------------------------------------
#### <a name="goals">Goals

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
- Include a tutorial

[Index](#top)


#### <a name="structure">Structure

- language(s): C++, python front end (c-extension? cython? need to choose how
  exactly we want to do this) [Sam]
- third party libraries: deal-II for FEM; PyNE for data?
- testing framework: [Google Test](https://github.com/google/googletest) 
- integration framework: [Travis CI](https://travis-ci.org/)
- test problem library (benchmark suite)
- build system: [Make](https://www.gnu.org/software/make/)
- documentation: [Doxygen](http://www.stack.nl/~dimitri/doxygen/)
- style guide: TBD [Josh]
- Collaboration / interaction plan
- License [slaybaugh]
- Anything else?

[Index](#top)


#### <a name="todo">To Do

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
- other transport models


[Index](#top)


