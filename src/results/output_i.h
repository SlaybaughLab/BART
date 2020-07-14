#ifndef BART_SRC_RESULTS_OUTPUT_I_H_
#define BART_SRC_RESULTS_OUTPUT_I_H_

#include <ostream>
#include <string>
#include <vector>

namespace bart {

namespace system {
struct System;
} // namespace system

namespace results {

class OutputI {
 public:
  virtual ~OutputI() = default;
  virtual void AddData(system::System &to_output) = 0;
  virtual void WriteData(std::ostream &output_stream) const = 0;
  virtual void WriteMasterFile(std::ostream &output_stream,
                               std::vector<std::string> filenames) const = 0;
  virtual void WriteVector(std::ostream &output_stream,
                           const std::vector<double> to_write) const = 0;
  virtual void WriteVector(std::ostream &output_stream,
                           const std::vector<double> to_write,
                           const std::vector<std::string> headers) const = 0;
};

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_I_H_
