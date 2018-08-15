#include <unistd.h>
#include <getopt.h>

#include <deal.II/base/mpi.h>
#include "gtest/gtest.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/bart_test_helper.h"

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
  return RUN_ALL_TESTS();
}

