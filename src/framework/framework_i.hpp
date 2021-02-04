#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_I_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_I_HPP_

#include <ostream>

#include "system/system.hpp"

namespace bart::framework {

class FrameworkI {
 public:
  virtual ~FrameworkI() = default;
  virtual void SolveSystem() = 0;
  virtual system::System* system() const = 0;
  virtual void OutputResults(std::ostream& output_stream) = 0;
  virtual void OutputMasterFile(std::ostream& output_stream,
                                const std::vector<std::string>& filenames,
                                const int process_id) = 0;
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_I_HPP_
