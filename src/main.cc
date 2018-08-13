#ifdef TEST
#include <unistd.h>
#include <getopt.h>

#include <deal.II/base/mpi.h>
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test_helpers/bart_test_helper.h"
#else
#include "common/problem_definition.h"
#include "common/bart_driver.h"
#endif



int main(int argc, char* argv[]) {
#ifdef TEST
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
#else
  try {
    if (argc!=2) {
      std::cerr << "Call the program as mpirun -np num_proc xtrans input_file_name" << std::endl;
      return 1;
    }
    ParameterHandler prm;
    bparams::DeclareParameters (prm);
    prm.read_input(argv[1]);
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    int dim = prm.get_integer ("problem dimension");
    switch (dim) {
      case 1:{
        BARTDriver<1> drive1d (prm);
        drive1d.DriveBART();
        break;
      }
      case 2:{
        BARTDriver<2> drive2d (prm);
        drive2d.DriveBART ();
        break;
      }
      case 3:{
        BARTDriver<3> drive3d (prm);
        drive3d.DriveBART ();
        break;
      }
      default:
        break;
    }
  }
  catch (std::exception &exc) {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
#endif
}
