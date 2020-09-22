#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string/to_string_converter.h"
#include "utility/colors.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

enum class StringColorPairToStringOutputTerm {
  kColorCode = 0, kString = 1, kColorReset = 2
};

class StringColorPairToString :
 public ToStringConverter<std::pair<std::string, utility::Color>,
                          StringColorPairToStringOutputTerm> {
 public:
  using OutputTerm = StringColorPairToStringOutputTerm;
  using Color = utility::Color;

  StringColorPairToString();
  virtual ~StringColorPairToString() = default;
  std::string Convert(const std::pair<std::string, Color> &) const override;
};

} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_STRING_COLOR_PAIR_TO_STRING_H_
