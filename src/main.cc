#ifdef TEST

#include <unistd.h>
#include <getopt.h>

#include "gtest/gtest.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/bart_test_helper.h"

#endif

//#include "aqdata/aq_base.h"

#ifdef TEST
int main(int argc, char* argv[]) {
  // Parse optional arguments
  ::testing::InitGoogleMock(&argc, argv);
  
  int option_index = 0;
  
  const struct option longopts[] =
  {
    {"report", no_argument, nullptr, 'r'}
  };
  
  int c;
  while ((c = getopt_long (argc, argv, "rd:", longopts, &option_index)) != -1) {
    switch(c) {
      case 'r':
        btest::GlobalBARTTestHelper().SetReport(true);
        break;
      case 'd':
        btest::GlobalBARTTestHelper().SetGoldFilesDirectory(optarg);
        break;
      default:
        break;
    }
  }
  //dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  // Testing
  //::testing::TestEventListeners& listeners =
  //::testing::UnitTest::GetInstance()->listeners();
  //if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) != 0)
  //if (world.rank() != 0)
  //delete listeners.Release(listeners.default_result_printer());
  return RUN_ALL_TESTS();
#else
int main() {
  return 0;
#endif
}
