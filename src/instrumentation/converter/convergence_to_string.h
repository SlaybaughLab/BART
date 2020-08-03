#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_

#include <map>
#include <string>
#include <variant>
#include <vector>

#include "instrumentation/converter/to_string_converter.h"
#include "convergence/status.h"

namespace bart {

namespace instrumentation {

namespace converter {

enum ConvergenceToStringOutputTerm {
  kIterationNum, kIterationMax, kDelta, kIndex
};

class ConvergenceToString :
    public ToStringConverter<convergence::Status, ConvergenceToStringOutputTerm> {
 public:
  using OutputTerm = ConvergenceToStringOutputTerm;
  using OutputTermToStringMap = std::map<OutputTerm, std::string>;

  ConvergenceToString()
      : ToStringConverter<convergence::Status, ConvergenceToStringOutputTerm>(
      "Iteration: ${ITERATION_NUM}/${ITERATION_MAX}, delta: ${DELTA}, index: ${INDEX}\n",
      {{kIterationNum, "${ITERATION_NUM}"},
       {kIterationMax, "${ITERATION_MAX}"},
       {kDelta, "${DELTA}"},
       {kIndex, "${INDEX}"}}) {}
  virtual ~ConvergenceToString() = default;
  std::string Convert(const convergence::Status& to_convert) const override;

  std::string SetOutputFormat(
      const std::vector<std::variant<OutputTerm, std::string>>);

  std::string null_character() const { return null_character_; };

 private:
  std::string null_character_{"∅"};
};



} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
