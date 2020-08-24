#include "utility/runtime/runtime_helper.h"

#include <getopt.h>
#include <unistd.h>
#include <iostream>

#include <deal.II/base/mpi.h>

namespace bart {

namespace utility {

namespace runtime {

namespace {
namespace MPI = dealii::Utilities::MPI;
}

RuntimeHelper::RuntimeHelper(std::string version)
    : n_mpi_processes_(MPI::n_mpi_processes(MPI_COMM_WORLD)),
      this_mpi_process_(MPI::this_mpi_process(MPI_COMM_WORLD)),
      version_(version) {}

void RuntimeHelper::ParseArguments(int argc, char **argv) {
  //int option_index = 0, c = 0;
  int c = 0;
  optind = 0;
  const struct option longopts[] = {
      {"pause",     no_argument,       NULL, 'p'},
      {NULL,        0,                 NULL,   0}
  };

  while ((c = getopt_long(argc, argv, "p", longopts, NULL)) != -1) {
    switch(c) {
      case 'p': {
        do_pause_ = true;
        break;
      }
    default:
        break;
    }
  }

  if (optind >= argc) {
    throw(std::runtime_error("No filename provided."));
  } else {
    filename_ = argv[optind];
  }
}


std::string RuntimeHelper::ProgramHeader() const {
  return std::string{"    __               __ \n"
                     "   / /_  ____  _____/ /_\n"
                     "  / __ \\/ __ `/ ___/ __/\n"
                     " / /_/ / /_/ / /  / /_  \n"
                     "/_.___/\\__,_/_/   \\__/  \n"
                     "BAY AREA RADIATION TRANSPORT\n"
                     "Developed at the University of California, Berkeley\n"
                     "version: " + version_ + "\n"};
}

} // namespace runtime

} // namespace utility

} // namespace bart
