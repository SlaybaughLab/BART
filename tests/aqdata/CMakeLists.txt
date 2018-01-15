ADD_EXECUTABLE(dummy_test_exe dummy_test.cc)
DEAL_II_SETUP_TARGET(dummy_test_exe)

SET(TEST_LIBRARIES _aqdata)
DEAL_II_PICKUP_TESTS()
