#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string_converter.h"

namespace bart {

namespace instrumentation {

namespace converter {

enum IntDoublePairToStringOutputTerm {kIndex, kValue};

class IntDoublePairToString
    : public ToStringConverter<std::pair<int, double>,
                               IntDoublePairToStringOutputTerm> {
 public:
  using OutputTerm = IntDoublePairToStringOutputTerm;

  IntDoublePairToString();
  virtual ~IntDoublePairToString() = default;

  std::string Convert(const std::pair<int, double> &input) const override;

  int precision() const {return precision_;}
 protected:
  int precision_ = 2;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_
