# Add a new executable for every single new test
ADD_EXECUTABLE(mesh_01_exe mesh_01.cc)
DEAL_II_SETUP_TARGET(mesh_01_exe)

# Set the libraries we wanna test
SET(TEST_LIBRARIES _mesh)

# Link the libraries
TARGET_LINK_LIBRARIES(mesh_01_exe ${TEST_LIBRARIES})
DEAL_II_PICKUP_TESTS()
