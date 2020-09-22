#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_

#include "instrumentation/converter/to_string/to_string_converter.h"

#include <complex>
#include <vector>

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

enum class IntVectorComplexPairToStringOutputTerm {
  kInt = 0, kVector = 1
};

enum class IntVectorComplexPairToStringVectorOutputTerm {
  kIndex = 0, kReal = 1, kImag = 2
};

class IntVectorComplexPairToString
 : public ToStringConverter<std::pair<int, std::vector<std::complex<double>>>,
                            IntVectorComplexPairToStringOutputTerm> {
 public:
  using InputType = std::pair<int, std::vector<std::complex<double>>>;
  using OutputTerm = IntVectorComplexPairToStringOutputTerm;
  using VectorTerm = IntVectorComplexPairToStringVectorOutputTerm;

  IntVectorComplexPairToString();
  std::string Convert(const InputType &input) const override;

  std::string SetVectorOutputFormat(
      const std::vector<std::variant<VectorTerm, std::string>>);

  IntVectorComplexPairToString& set_vector_delimiter(const std::string to_set) {
    vector_delimiter_ = to_set;
    return *this;
  }
  IntVectorComplexPairToString& set_precision(const int to_set) {
    precision_ = to_set;
    return *this; }

  int precision() const { return precision_; }
  std::string vector_delimiter() const { return vector_delimiter_; }
  std::string vector_output_format() const { return vector_output_format_; }
  std::map<VectorTerm, std::string> vector_term_to_string_map() const {
    return vector_term_to_string_map_; }
 private:
  int precision_{2};
  std::string vector_delimiter_{", "};
  std::string vector_output_format_{""};
  std::map<VectorTerm, std::string> vector_term_to_string_map_;
};

} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart


#endif //BART_SRC_INSTRUMENTATION_CONVERTER_INT_VECTOR_COMPLEX_PAIR_TO_STRING_H_
