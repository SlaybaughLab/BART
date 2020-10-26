#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string/to_string_converter.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

enum class IntDoublePairToStringOutputTerm {kIndex = 0, kValue = 1};

class IntDoublePairToString
    : public ToStringConverter<std::pair<int, double>,
                               IntDoublePairToStringOutputTerm> {
 public:
  using OutputTerm = IntDoublePairToStringOutputTerm;

  IntDoublePairToString();
  virtual ~IntDoublePairToString() = default;

  std::string Convert(const std::pair<int, double> &input) const override;

  IntDoublePairToString& set_precision(const int to_set) {
    precision_ = to_set;
    return *this; }
  int precision() const {return precision_;}
 protected:
  int precision_ = 2;
  static bool is_registered_;
};

} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_INT_DOUBLE_PAIR_TO_STRING_H_
