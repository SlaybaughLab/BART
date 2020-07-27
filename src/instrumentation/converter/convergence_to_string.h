#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_

#include <map>
#include <string>

#include "instrumentation/converter/converter_i.h"
#include "convergence/status.h"

namespace bart {

namespace instrumentation {

namespace converter {

class ConvergenceToString : public ConverterI<convergence::Status, std::string> {
 public:
  enum OutputTerm {
    kIterationNum, kIterationMax, kDelta, kIndex
  };
  using OutputTermToStringMap = std::map<OutputTerm, std::string>;

  virtual ~ConvergenceToString() = default;
  std::string Convert(const convergence::Status& to_convert) const override {
    return std::string{};
  };

  std::string output_format() const { return output_format_; }
  OutputTermToStringMap output_term_to_string_map() const {
    return output_term_to_string_map_; }


 private:
  std::string output_format_{"Iteration: ITERATION_NUM/ITERATION_MAX, delta: DELTA, index: INDEX"};
  OutputTermToStringMap output_term_to_string_map_{
      {kIterationNum, "ITERATION_NUM"},
      {kIterationMax, "ITERATION_MAX"},
      {kDelta, "DELTA"},
      {kIndex, "INDEX"}
  };
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
