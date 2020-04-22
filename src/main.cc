#include <memory>
#include <iostream>
#include <fstream>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/mpi.h>

#include "framework/builder/framework_builder.h"
#include "utility/reporter/mpi.h"
#include "problem/parameters_dealii_handler.h"

int main(int argc, char* argv[]) {
  try {
    if (argc != 2) {
      std::cerr
          << "Call the program as mpirun -np num_proc bart input_file_name"
          << std::endl;
      return 1;
    }

    std::cout << "BAY AREA RADIATION TRANSPORT\n"
              << "Developed at the University of California, Berkeley"
              << std::endl;

    bart::problem::ParametersDealiiHandler prm;
    dealii::ParameterHandler d2_prm;
    const std::string filename{argv[1]};

    prm.SetUp(d2_prm);
    d2_prm.parse_input(filename, "");
    prm.Parse(d2_prm);

    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    double k_eff_final;

    const int n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const int process_id = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    std::ofstream output_stream;
    const std::string output_filename_base{prm.OutputFilenameBase()};
    if (n_processes > 1) {
      const std::string full_filename = output_filename_base + dealii::Utilities::int_to_string(process_id, 4);
      output_stream.open((full_filename + ".vtu").c_str());
    } else {
      output_stream.open((output_filename_base + ".vtu").c_str());
    }

    auto reporter_ptr = std::make_shared<bart::utility::reporter::Mpi>(
        std::make_unique<dealii::ConditionalOStream>(std::cout, process_id == 0));
    switch(prm.SpatialDimension()) {
      case 1: {
        bart::framework::builder::FrameworkBuilder<1> builder(reporter_ptr);
        auto framework_ptr = builder.BuildFramework("main", prm);
        framework_ptr->SolveSystem();
        framework_ptr->OutputResults(output_stream);
        k_eff_final = framework_ptr->system()->k_effective.value_or(0);
        break;
      }
      case 2: {
        bart::framework::builder::FrameworkBuilder<2> builder(reporter_ptr);
        auto framework_ptr = builder.BuildFramework("main", prm);
        framework_ptr->SolveSystem();
        framework_ptr->OutputResults(output_stream);
        k_eff_final = framework_ptr->system()->k_effective.value_or(0);
        break;
      }
      case 3: {
        bart::framework::builder::FrameworkBuilder<3> builder(reporter_ptr);
        auto framework_ptr = builder.BuildFramework("main", prm);
        framework_ptr->SolveSystem();
        framework_ptr->OutputResults(output_stream);
        k_eff_final = framework_ptr->system()->k_effective.value_or(0);
        break;
      }
    }

    std::cout << "Final k_effective: " << k_eff_final << std::endl;

  } catch (std::exception &exc) {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
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
