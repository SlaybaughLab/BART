#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string_converter.h"

#include <complex>
#include <vector>

namespace bart {

namespace instrumentation {

namespace converter {

enum class IntVectorComplexPairToStringOutputTerm {
  kInt = 0, kVectorIndex = 1, kReal = 2, kImag = 3
};

class IntVectorComplexPairToString
 : public ToStringConverter<std::pair<int, std::vector<std::complex<double>>>,
                            IntVectorComplexPairToStringOutputTerm> {
 public:
  using InputType = std::pair<int, std::vector<std::complex<double>>>;
  using OutputTerm = IntVectorComplexPairToStringOutputTerm;

  IntVectorComplexPairToString();
  std::string Convert(const std::pair<int,
                                      std::vector<std::complex<double>>> &input) const override;

  IntVectorComplexPairToString& set_precision(const int to_set) {
    precision_ = to_set;
    return *this; }
  int precision() const { return precision_; }
 private:
  int precision_ = 2;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart


#endif //BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_
