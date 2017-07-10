## <a name="top">Meeting Notes

* [Goals](#goals)
* [Structure](#structure)
* [To Do](#todo)


-----------------------------------------------------------------
#### <a name="goals">Goals

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
- be able to handle small and medium problems (up to ~100M degrees of freedom)

[Index](#top)


#### <a name="structure">Structure

- language(s): C++, Python 3 front end (c-extension? cython? need to choose how
  exactly we want to do this) [Sam]
- third party libraries: 
  - [deal-II](http://www.dealii.org/) for FEM
  - [PyNE](https://github.com/pyne/pyne) modules pulled in for nuclear data
- testing framework: [Google Test](https://github.com/google/googletest) 
- integration framework: [Travis CI](https://travis-ci.org/)
- test problem library (benchmark suite)
- build system: [Make](https://www.gnu.org/software/make/)
- documentation: [Doxygen](http://www.stack.nl/~dimitri/doxygen/)
- style guide: TBD [Josh]
- Collaboration / interaction plan
- [MIT License](https://github.com/SlaybaughLab/BART/blob/master/LICENSE)
- Anything else?

[Index](#top)


#### <a name="todo">To Do

**near-term** (next 6 months)

Goal: structure in place; two methods available

Modules:
- implement with testing, CI, documentation, etc. 
- Cartesian spatial solver (2D and/or 3D)
- access to nuclear data
- source specification
  - input 
  - internal representation
- geometry specification
- mesh: good idea to be able to take [Exodous II](https://cubit.sandia.gov/public/13.2/help_manual/WebHelp/finite_element_model/exodus/exodus2_file_specification.htm) format
- boundary conditions (vacuum and incident)
- material assignment
- multiple energy groups

2017-07-06 Plan of action for code integration
- merge [DG-EP repo](https://github.com/weixiong-zheng-berkeley/DG-EP) with
  this repo {done}
- implement test framework
- hook into CI framework
- update documentation to doxygen


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


