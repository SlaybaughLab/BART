##
#  CMake project for DG-EP demonstration based on deal.II
##

# Set the name of the project and target:
SET(TARGET xtrans)

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
#  SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
# SET(TARGET_SRC
#  ${TARGET}.cc
#  )
FILE(GLOB TARGET_SRCS "src/*/*.*" "src/main.cc")
# FILE(GLOB TARGET_INCS "src/*/*.h")
# SET(TARGET_SRC ${TARGET_SRCS} ${TARGET_INCS})
SET(TARGET_SRC ${TARGET_SRCS})
# SET(TARGET_INC ${TARGET_INC})

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.4 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Are all dependencies fulfilled?
#
IF(NOT (DEAL_II_WITH_PETSC OR DEAL_II_WITH_TRILINOS) OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_P4EST = ON
or
    DEAL_II_WITH_TRILINOS = ON
    DEAL_II_WITH_P4EST = ON
One or both of these combinations are OFF in your installation but at least one is required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2 *.vtu *.pvtu)
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

# separating source files (excluding main.cc) to make a library
ADD_LIBRARY (xtrans_lib src/*/*.cc)
ADD_EXECUTABLE (xtrans src/main.cc)

# enable the testing
ENABLE_TESTING ()
FIND_PACKAGE (gtest REQUIRED)
FIND_PACKAGE (gmock REQUIRED)

FILE (GLOB TESTING_SRCS "src/*/tests/*.cc")
SET (TESTING_SRC_FILES ${TESTING_SRCS})
ADD_EXECUTABLE (xtrans_test, ${TESTING_SRC_FILES})
ADD_TEST (NAME unit-tests COMMAND xtrans_test)

# link with gtest gmock
TARGET_LINK_LIBRARIES (xtrans_test gtest gtest_main gmock gmock_main)

# link with xtrans library
TARGET_LINK_LIBRARIES (xtrans_test xtrans_lib)
