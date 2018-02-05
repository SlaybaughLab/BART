# Add a new executable for every single new test
ADD_EXECUTABLE(mpi_dummy_exe mpi_dummy.cc)
DEAL_II_SETUP_TARGET(mpi_dummy_exe)

# Set the libraries we wanna test
SET(TEST_LIBRARIES _mesh)

# Link the libraries
TARGET_LINK_LIBRARIES(mpi_dummy_exe ${TEST_LIBRARIES})
DEAL_II_PICKUP_TESTS()
