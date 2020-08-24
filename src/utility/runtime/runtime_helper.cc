#include "utility/runtime/runtime_helper.h"

#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>

namespace bart {

namespace utility {

namespace runtime {

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
                     "   / /_  ____ ______/ /_\n"
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
