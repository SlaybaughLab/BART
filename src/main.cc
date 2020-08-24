#include <memory>
#include <iostream>
#include <fstream>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/mpi.h>

#include "framework/builder/framework_builder.h"
#include "problem/parameters_dealii_handler.h"
#include "utility/runtime/runtime_helper.h"

int main(int argc, char* argv[]) {
  try {
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    bart::utility::runtime::RuntimeHelper runtime_helper("0.2.2");
    std::cout << runtime_helper.ProgramHeader() << std::endl;
    runtime_helper.ParseArguments(argc, argv);
    const std::string filename{runtime_helper.filename()};
    const int n_processes = runtime_helper.n_mpi_processes();
    const int process_id =  runtime_helper.this_mpi_process();

    bart::problem::ParametersDealiiHandler prm;
    dealii::ParameterHandler d2_prm;

    prm.SetUp(d2_prm);
    d2_prm.parse_input(filename, "");
    prm.Parse(d2_prm);

    double k_eff_final;

    // Open file for output, if there are multiple processes they will end with
    // a number indicating the process number.
    std::ofstream output_stream, master_output_stream;
    const std::string output_filename_base{prm.OutputFilenameBase()};
    if (n_processes > 1) {
      const std::string full_filename = output_filename_base + dealii::Utilities::int_to_string(process_id, 4);
      output_stream.open((full_filename + ".vtu").c_str());
      master_output_stream.open(output_filename_base + ".pvtu");
    } else {
      output_stream.open((output_filename_base + ".vtu").c_str());
    }
    std::vector<std::string> filenames;
    // Write master pvtu record if multiple files are required
    if ((n_processes > 1) && (process_id == 0)) {
      for (int process = 0; process < n_processes; ++process) {
        const std::string full_filename =
            output_filename_base + dealii::Utilities::int_to_string(process, 4);
        filenames.push_back(full_filename + ".vtu");
      }
    }

    // Framework pointer
    std::unique_ptr<bart::framework::FrameworkI> framework_ptr;

    switch(prm.SpatialDimension()) {
      case 1: {
        bart::framework::builder::FrameworkBuilder<1> builder;
        framework_ptr = builder.BuildFramework("main", prm);
        break;
      }
      case 2: {
        bart::framework::builder::FrameworkBuilder<2> builder;
        framework_ptr = builder.BuildFramework("main", prm);
        break;
      }
      case 3: {
        bart::framework::builder::FrameworkBuilder<3> builder;
        framework_ptr = builder.BuildFramework("main", prm);
        break;
      }
    }

    // Pause if needed before solve
    if (runtime_helper.do_pause()) {
      std::cout << "Press <Enter> to begin solve...";
      std::cin.ignore();
    }

    framework_ptr->SolveSystem();
    framework_ptr->OutputResults(output_stream);
    if (n_processes > 1)
      framework_ptr->OutputMasterFile(master_output_stream, filenames, process_id);
    k_eff_final = framework_ptr->system()->k_effective.value_or(0);

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
