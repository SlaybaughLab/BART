#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <vector>


#include <deal.II/base/mpi.h>
#include "gtest/gtest.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/bart_test_helper.h"

int main(int argc, char** argv) {
  // Parse optional arguments

  int option_index = 0, c = 0, listener = -1;
  bool use_mpi = false;
  std::string filter = "";

  const struct option longopts[] =
  {
    {"filter",    no_argument,       NULL, 'f'},
    {"report",    no_argument,       NULL, 'r'},
    {"mpi",       no_argument,       NULL, 'm'},
    {"listeners", required_argument, NULL, 'l'},
    {NULL,        0,                 NULL,   0}
  };

  while ((c = getopt_long (argc, argv, "rmd:f:l:", longopts, &option_index)) != -1) {
    switch(c) {
      case 'r': {
        bart::test_helpers::GlobalBARTTestHelper().SetReport(true);
        break;
      }
      case 'd': {
        bart::test_helpers::GlobalBARTTestHelper().SetGoldFilesDirectory(optarg);
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
      case 'l': {
        listener = atoi(optarg);
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
    new_argv.push_back("--gtest_filter=-*MPIOnly*");
  }

  new_argv.push_back(nullptr);
  argv = const_cast<char**>(new_argv.data());
  argc += 1;

  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  ::testing::InitGoogleMock(&argc, argv);

  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (listener != -1) {
    if (world_rank != listener) {
      delete listeners.Release(listeners.default_result_printer());
    }
  }
  
  // Re-seed random number generator for random number tests
  std::srand(time(NULL));
  return RUN_ALL_TESTS();
}

