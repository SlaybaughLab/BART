##
#  Demo for separating CMakeLists.txt
##

# Set the name of the project and target:
# SET(TARGET xtrans)

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
# FILE(GLOB TARGET_SRC "src/*/*.cc")
# FILE(GLOB TARGET_INC "src/*/*.h")
# SET(TARGET_SRC ${TARGET_SRC})
# SET(TARGET_INC ${TARGET_INC})

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(xtrans)

### DEPENDENCIES #####################################################
# Check that DEAL II is installed
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

### GTEST ############################################################

find_package(GTest REQUIRED)
set(GMOCK_BOTH_LIBRARIES "/usr/lib/libgmock.so;/usr/lib/libgmock_main.so")
include_directories(${GTEST_INCLUDE_DIRS})

######################################################################

DEAL_II_INITIALIZE_CACHED_VARIABLES()
SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd
  *.d2 *.vtu *.pvtu)
  
### EXECUTABLES AND SOURCES ##########################################

# Recursively find source files

file(GLOB_RECURSE sources "src/*.cpp" "src/*.cc")
set(testing_sources ${sources})
list(FILTER sources EXCLUDE REGEX ".(/tests/).")

# Add xtrans executables
ADD_EXECUTABLE(xtrans ${sources})
ADD_EXECUTABLE(xtrans_test ${testing_sources})

# Add testing definition and library to xtrans_test
target_compile_definitions(xtrans_test PUBLIC -DTEST)
target_link_libraries(xtrans_test ${GTEST_BOTH_LIBRARIES})

DEAL_II_SETUP_TARGET(xtrans)
DEAL_II_SETUP_TARGET(xtrans_test)

### COPY TEST FILES ##################################################

file(GLOB_RECURSE test_files "src/*.golden")
foreach(test_file ${test_files})
  MESSAGE(${test_file})
  get_filename_component(file_name ${test_file} NAME)
  add_custom_command(
    TARGET xtrans_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${test_file}
            ${CMAKE_CURRENT_BINARY_DIR}/test_data/${file_name})
endforeach()  
