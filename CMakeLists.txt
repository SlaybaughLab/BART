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

### GTEST/GMOCK  ######################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(GTEST_ROOT "/usr")
find_package(GTest REQUIRED 1.8.0)
find_package(GMock REQUIRED)
#set(GTEST_INCLUDE_DIRS "/usr/include")
#set(GTEST_BOTH_LIBRARIES "/usr/lib/libgtest_main.so;/usr/lib/libgtest.so")
#set(GMOCK_BOTH_LIBRARIES "/usr/lib/libgmock.a;/usr/lib/libgmock_main.a")

### COMPILE AND LINKER FLAGS #########################################

SET(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -g -O0 --coverage")
SET(DEAL_II_LINKER_FLAGS_DEBUG "${DEAL_II_LINKER_FLAGS_DEBUG} --coverage")

######################################################################

DEAL_II_INITIALIZE_CACHED_VARIABLES()
SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd
  *.d2 *.vtu *.pvtu)
ENABLE_TESTING()

### EXECUTABLES AND SOURCES ##########################################

# Recursively find source files

file(GLOB_RECURSE sources "src/*.cpp" "src/*.cc")
set(testing_sources ${sources})
list(FILTER sources EXCLUDE REGEX ".(/tests/).")
list(FILTER sources EXCLUDE REGEX ".(/test_helpers/).")

# Include directories
include_directories(${GTEST_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/inc)

# Add xtrans executables
ADD_EXECUTABLE(xtrans ${sources})
ADD_EXECUTABLE(xtrans_test ${testing_sources})

# Add testing definition and library to xtrans_test
target_compile_definitions(xtrans_test PUBLIC -DTEST)

target_link_libraries(xtrans_test ${GTEST_BOTH_LIBRARIES}
  ${GMOCK_BOTH_LIBRARIES})

# Add subdirectory builds to build libraries for CTEST
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/aqdata)
#ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/common)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/mesh)

SET(TARGET "xtrans")
SET(TARGET_SRC ${CMAKE_SOURCE_DIR}/src/main.cc)

# Add unit testing subdirectories
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/tests/aqdata)
#ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/tests/common)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/tests/mesh)

DEAL_II_SETUP_TARGET(xtrans)
DEAL_II_SETUP_TARGET(xtrans_test)

### TEST FILES ##################################################
# Create copies of .gold files in the BART/src directory for gtest

add_custom_target(copy_gtest_gold_files ALL)

add_custom_command(TARGET copy_gtest_gold_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory test_data)

file(GLOB_RECURSE gtest_files "src/*.gold")
foreach(test_file ${gtest_files})
  get_filename_component(file_name ${test_file} NAME)
  add_custom_command(
    TARGET copy_gtest_gold_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${test_file}
            ${CMAKE_CURRENT_BINARY_DIR}/test_data/${file_name})
endforeach()
add_dependencies(xtrans_test copy_gtest_gold_files)

# Creates copies of .output files in the BART/test directory for CTEST

add_custom_target(copy_ctest_output_files ALL)

file(GLOB_RECURSE ctest_files "tests/*.output")
foreach(test_file ${ctest_files})
  string(REGEX REPLACE ".*BART/tests" "${CMAKE_CURRENT_BINARY_DIR}/tests"
    file_name ${test_file})
  add_custom_command(
    TARGET copy_ctest_output_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${test_file} ${file_name})
endforeach()
add_dependencies(xtrans_test copy_ctest_output_files)
