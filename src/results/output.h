#ifndef BART_SRC_RESULTS_OUTPUT_H_
#define BART_SRC_RESULTS_OUTPUT_H_

#include "results/output_i.h"

namespace bart {

namespace results {

class Output : public OutputI {
 public:
  virtual ~Output() = default;
  void WriteVector(std::ostream &output_stream,
                           const std::vector<double> to_write) const override;

  void WriteVector(std::ostream &output_stream,
                   std::vector<double> to_write,
                   std::vector<std::string> headers) const override;
};

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_H_
