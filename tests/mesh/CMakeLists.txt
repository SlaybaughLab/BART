# Add a new executable for every single new test
ADD_DEFINITIONS(-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(mpi_dummy_exe mpi_dummy.cc)
DEAL_II_SETUP_TARGET(mpi_dummy_exe)

ADD_EXECUTABLE(mesh_generator_01_exe mesh_generator_01.cc)
DEAL_II_SETUP_TARGET(mesh_generator_01_exe)
# Set the libraries we wanna test
SET(TEST_LIBRARIES _mesh)

# Link the libraries
TARGET_LINK_LIBRARIES(mpi_dummy_exe ${TEST_LIBRARIES})
TARGET_LINK_LIBRARIES(mesh_generator_01_exe ${TEST_LIBRARIES})
DEAL_II_PICKUP_TESTS()
