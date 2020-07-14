#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_I_H_
#define BART_SRC_FRAMEWORK_FRAMEWORK_I_H_

#include <ostream>

#include "system/system.h"

namespace bart {

namespace framework {

class FrameworkI {
 public:
  virtual ~FrameworkI() = default;
  virtual void SolveSystem() = 0;
  virtual system::System* system() const = 0;
  virtual void OutputResults(std::ostream& output_stream) = 0;
  virtual void OutputMasterFile(std::ostream& output_stream,
                                const std::vector<std::string>& filenames,
                                const int process_id) = 0;
  virtual void OutputIterationError(std::ostream& output_stream) = 0;
};

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_I_H_
