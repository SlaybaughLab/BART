#include "../test_utilities.h"

void test ()
{
  dealii::deallog << "Process ID: "
                  << dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)
                  << std::endl;
}

int main (int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  testing::MPILogInit init;
  dealii::deallog.push ("MPI Demo");
  test ();
  dealii::deallog.pop ();
  return 0;
}
