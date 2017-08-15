## <a name="top">Meeting Notes

* [Goals](#goals)
* [Structure](#structure)
* [To Do](#todo)
* [To Do July](#july)

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
- covariant reference frame [Sam]
- include a tutorial
- be able to handle small and medium problems (up to ~100M degrees of freedom)

[Index](#top)


-----------------------------------------------------------------
#### <a name="structure">Structure

- language(s): C++, Python 3 front end (c-extension? cython? need to choose how
  exactly we want to do this) [Sam]
- third party libraries: 
  - [deal-II](http://www.dealii.org/) for FEM
  - [PyNE](https://github.com/pyne/pyne) modules pulled in for nuclear data
- build environment: [Docker](https://www.docker.com/)?
- testing framework: [Google Test](https://github.com/google/googletest) 
- integration framework: [GoCD](https://www.gocd.org/)?  ~~[Travis CI](https://travis-ci.org/)~~
- test problem library (benchmark suite)
- build system: [Cmake](https://cmake.org/)
- documentation: [Doxygen](http://www.stack.nl/~dimitri/doxygen/)
- style guide: [Google C++](https://google.github.io/styleguide/cppguide.html)
- Collaboration / interaction plan: 
  - each new feature will get its own branch until it is ready to be merged into
    master
  - code contributions should be submitted as a PR to the main repo, someone
    else should review the code and merge it
  - bug / feature / etc. requests should start as issues
  - [rebase](https://git-scm.com/docs/git-rebase): use this strategy for
    branches; 
    [here](https://www.atlassian.com/git/tutorials/merging-vs-rebasing) is some info
    about how and why. Note: use with caution.
  - we will use tags for major releases
- [MIT License](https://github.com/SlaybaughLab/BART/blob/master/LICENSE)
- Anything else?

[Index](#top)


-----------------------------------------------------------------
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
- mesh: good idea to be able to take [Exodous II](https://cubit.sandia.gov/public/13.2/help_manual/WebHelp/finite_element_model/exodus/exodus2_file_specification.htm) format; this works with Cubit/[Trelis](http://www.csimsoft.com/trelis), and we have two Trelis licenses.
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



-----------------------------------------------------------------
#### <a name="july">To Do for the month of July
**2017-07-13**

Notes about structure:
- Header and source code files need to go in the same folder
- Each class gets its own folder
- Derived classes can go in the same folder as the base class
- Tests for each piece of code should be in a `test` folder in the directory
  where the code it is testing lives
- We will have these folders
  - `src`, which will contain all code and integration tests
  - `doc`, which will contain a `development` folder and a `code_doc` foder
- We will use linux and osx for testing on Travis CI [info
  here](https://docs.travis-ci.com/user/multi-os/)


**2017-07-06** 
- ~~merge [DG-EP repo](https://github.com/weixiong-zheng-berkeley/DG-EP) with~~
  ~~this repo {done}~~
- implement test framework
- hook into CI framework
- update documentation to doxygen


[Index](#top)
