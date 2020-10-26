#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_PAIR_INCREMENTER_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_PAIR_INCREMENTER_H_

#include "instrumentation/converter/converter_i.h"
#include "instrumentation/converter/factory.hpp"

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType>
class PairIncrementer : public ConverterI<InputType, std::pair<int, InputType>> {
 public:
  std::pair<int, InputType> Convert(const InputType &input) const override;
 private:
  mutable int increment_{-1};
  static bool is_registered_;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_PAIR_INCREMENTER_H_
