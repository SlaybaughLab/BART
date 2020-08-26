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
    : version_(version) {}

void RuntimeHelper::ParseArguments(int argc, char **argv) {
  //int option_index = 0, c = 0;
  int c = 0;
  optind = 0;
  const struct option longopts[] = {
      {"help",      no_argument,       NULL, 'h'},
      {"pause",     no_argument,       NULL, 'p'},
      {NULL,        0,                 NULL,   0}
  };

  while ((c = getopt_long(argc, argv, "hp", longopts, NULL)) != -1) {
    switch(c) {
      case 'p': {
        do_pause_ = true;
        break;
      }
      case '?': {
        std::cerr << "Unknown argument: " + std::to_string(optopt) + "\n\n";
        [[fallthrough]];
      }
      case 'h':
      default:
        show_help_ = true;
        break;
    }
  }

  if (optind >= argc && !show_help_) {
    throw(std::runtime_error("No filename provided. See bart --help for usage."));
  } else if(!show_help_) {
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
std::string RuntimeHelper::HelpMessage() const {
  return std::string{"usage: bart [OPTION] <filename>\n"
                     "  -h, --help\t\tDisplays this help text.\n"
                     "  -p, --pause\t\tPauses after building frameworks and "
                     "before commencing solve."};
}

} // namespace runtime

} // namespace utility

} // namespace bart
