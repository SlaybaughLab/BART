#include "utility/runtime/runtime_helper.h"

namespace bart {

namespace utility {

namespace runtime {

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
