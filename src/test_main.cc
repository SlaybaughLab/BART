#include <unistd.h>
#include <getopt.h>
#include <vector>

#include <deal.II/base/mpi.h>
#include "gtest/gtest.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/bart_test_helper.h"

int main(int argc, char** argv) {
  // Parse optional arguments

  int option_index = 0, c = 0;
  bool use_mpi = false;
  std::string filter = "";

  const struct option longopts[] =
  {
    {"filter", no_argument, NULL, 'f'},
    {"report", no_argument, NULL, 'r'},
    {"mpi",    no_argument, NULL, 'm'},
    {NULL,     0,           NULL,   0}
  };

  while ((c = getopt_long (argc, argv, "rmd:f:", longopts, &option_index)) != -1) {
    switch(c) {
      case 'r': {
        btest::GlobalBARTTestHelper().SetReport(true);
        break;
      }
      case 'd': {
        btest::GlobalBARTTestHelper().SetGoldFilesDirectory(optarg);
        break;
      }
      case 'm': {
        use_mpi = true;
        break;
      }
      case 'f': {
        filter = optarg;
        break;
      }
      default:
        break;
    }
  }

  std::vector<const char*> new_argv(argv, argv + argc);
  if (filter != "") {
    filter = "--gtest_filter=" + filter;
    new_argv.push_back(filter.c_str());
  } else if (use_mpi) {
      new_argv.push_back("--gtest_filter=*MPI*");
  } else {
    new_argv.push_back("--gtest_filter=-*MPI*");
  }

  new_argv.push_back(nullptr);
  argv = const_cast<char**>(new_argv.data());
  argc += 1;

  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  ::testing::InitGoogleMock(&argc, argv);
  
  return RUN_ALL_TESTS();
}

