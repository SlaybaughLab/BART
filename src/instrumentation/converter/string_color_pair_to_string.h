#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string_converter.h"
#include "utility/colors.h"

namespace bart {

namespace instrumentation {

namespace converter {

enum class StringColorPairToStringOutputTerm {
  kColorCode = 0, kString = 1, kColorReset = 2
};

class StringColorPairToString :
 public ToStringConverter<std::pair<std::string, utility::Color>,
                          StringColorPairToStringOutputTerm> {
 public:
  using OutputTerm = StringColorPairToStringOutputTerm;

  StringColorPairToString();
  virtual ~StringColorPairToString() = default;
  std::string Convert(const std::pair<std::string,
                                      utility::Color> &input) const override;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_
