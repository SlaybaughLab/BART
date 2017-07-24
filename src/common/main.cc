/* ---------------------------------------------------------------------
 *
 *
 * Author: Weixiong Zheng
 *
 *
 * ----------------------------------------------------------------------
 */
//#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

#include "problem_definition.h"
#include "bart_driver.h"

using namespace dealii;

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;
    
    if (argc!=2)
    {
      std::cerr << "Call the program as mpirun -np num_proc xtrans input_file_name" << std::endl;
      return 1;
    }
    ParameterHandler prm;
    ProblemDefinition::declare_parameters (prm);
    prm.read_input(argv[1]);
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    // ModelManager modeler (prm);
    dimension =
    modeler.build_and_run_model (prm);
  }
  catch (std::exception &exc)
  {
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
  catch (...)
  {
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
}
