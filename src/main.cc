#ifdef TEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <unistd.h>
#include <cstdlib>
#include <deal.II/base/mpi.h>

#include "test_helpers/bart_test_helper.h"

#endif

//#include "aqdata/aq_base.h"

#ifdef TEST
int main(int argc, char* argv[]) {
  // Parse optional arguments
  int c;
  while ((c = getopt (argc, argv, "r")) != -1)
    switch(c) {
      case 'r':
        btest::GlobalBartTestHelper().ReInit(true, "test_data/");
    }
  ::testing::InitGoogleMock(&argc, argv);
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  // Testing
  ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) != 0)
  //if (world.rank() != 0)
    delete listeners.Release(listeners.default_result_printer());
  return RUN_ALL_TESTS();
#else
int main() {
  return 0;
#endif
}
