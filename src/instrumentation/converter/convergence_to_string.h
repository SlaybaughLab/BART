#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_

#include <string>

#include "instrumentation/converter/converter_i.h"
#include "convergence/status.h"

namespace bart {

namespace instrumentation {

namespace converter {

class ConvergenceToString : public ConverterI<convergence::Status, std::string> {
 public:
  virtual ~ConvergenceToString() = default;
  std::string Convert(const convergence::Status& to_convert) const override {
    return std::string{};
  };
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_CONVERGENCE_TO_STRING_H_
