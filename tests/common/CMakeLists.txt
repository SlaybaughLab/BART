# Add a new executable for every single new test
ADD_EXECUTABLE(problem_definition_01_exe problem_definition_01.cc)
DEAL_II_SETUP_TARGET(problem_definition_01_exe)

ADD_EXECUTABLE(bart_builder_01_exe bart_builder_01.cc)
DEAL_II_SETUP_TARGET(bart_builder_01_exe)

# Set the libraries we wanna test
SET(TEST_LIBRARIES _common)

# Link the libraries
TARGET_LINK_LIBRARIES(problem_definition_01_exe ${TEST_LIBRARIES})
TARGET_LINK_LIBRARIES(bart_builder_01_exe ${TEST_LIBRARIES})

DEAL_II_PICKUP_TESTS()
