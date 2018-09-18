#
#  ____          _____ _______ 
# |  _ \   /\   |  __ \__   __|
# | |_) | /  \  | |__) | | |   
# |  _ < / /\ \ |  _  /  | |   
# | |_) / ____ \| | \ \  | |   
# |____/_/    \_\_|  \_\ |_|   
#                              
# Bay Area Radiation Transport
#
# CMakeLists derived from dealii example CMakeFiles
#
# Author: J.S. Rehak <jsrehak@berkeley.edu>
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(xtrans)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

### PROTOCOL BUFFERS #################################################

#include(FindProtobuf)
find_package(Protobuf REQUIRED)
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS
  "${CMAKE_SOURCE_DIR}/protobuf/material.proto")
get_filename_component(PROTO_INC_DIR ${PROTO_HDRS} DIRECTORY)
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
list(APPEND sources ${PROTO_SRCS} ${PROTO_HDRS})
set(testing_sources ${sources})
list(FILTER sources EXCLUDE REGEX ".*/tests/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_helpers/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_main.cc$")
list(FILTER testing_sources EXCLUDE REGEX ".*/main.cc$")

# Include directories
include_directories(${GTEST_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/inc
  ${PROTO_INC_DIR})

# Add xtrans executables
ADD_EXECUTABLE(xtrans ${sources})
ADD_EXECUTABLE(xtrans_test ${testing_sources})

# Add testing definition and library to xtrans_test
target_compile_definitions(xtrans_test PUBLIC -DTEST)


target_link_libraries(xtrans ${Protobuf_LIBRARIES})
target_link_libraries(xtrans_test ${GTEST_BOTH_LIBRARIES}
  ${GMOCK_BOTH_LIBRARIES} ${Protobuf_LIBRARIES})

# Add subdirectory builds to build libraries for CTEST
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/mesh)

# Add unit testing subdirectories
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

# Create copies of material protocol buffer data
add_custom_target(copy_protoc_files ALL)
file(COPY "./src/material/tests/data/" DESTINATION "./test_data/material")
add_dependencies(xtrans_test copy_protoc_files)
