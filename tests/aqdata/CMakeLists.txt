# Add a new executable for every single new test
ADD_EXECUTABLE(dummy_test_exe dummy_test.cc)
DEAL_II_SETUP_TARGET(dummy_test_exe)

ADD_EXECUTABLE(aq_base_01_exe aq_base_01.cc)
DEAL_II_SETUP_TARGET(aq_base_01_exe)

ADD_EXECUTABLE(aq_base_1d_exe aq_base_1d.cc)
DEAL_II_SETUP_TARGET(aq_base_1d_exe)

ADD_EXECUTABLE(lsgc_01_exe lsgc_01.cc)
DEAL_II_SETUP_TARGET(lsgc_01_exe)

ADD_EXECUTABLE(lsgc_non_ep_exe lsgc_non_ep.cc)
DEAL_II_SETUP_TARGET(lsgc_non_ep_exe)

# Set the libraries we wanna test
SET(TEST_LIBRARIES _aqdata)

# Link the libraries
TARGET_LINK_LIBRARIES(aq_base_01_exe ${TEST_LIBRARIES})
TARGET_LINK_LIBRARIES(aq_base_1d_exe ${TEST_LIBRARIES})
TARGET_LINK_LIBRARIES(lsgc_01_exe ${TEST_LIBRARIES})
TARGET_LINK_LIBRARIES(lsgc_non_ep_exe ${TEST_LIBRARIES})
DEAL_II_PICKUP_TESTS()
