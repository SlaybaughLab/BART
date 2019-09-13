#ifndef BART_SRC_RESULTS_OUTPUT_I_H_
#define BART_SRC_RESULTS_OUTPUT_I_H_

namespace bart {

namespace system {
struct System;
} // namespace system

namespace results {

class OutputI {
 public:
  ~OutputI() = default;
  virtual void Output(system::System &to_output) const = 0;
};

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_I_H_
