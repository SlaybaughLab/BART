#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

#include "instrumentation/converter/converter_i.h"

namespace bart {

namespace instrumentation {

namespace converter {

class DealiiToComplexVector :
    public ConverterI<dealii::Vector<double>,
                      std::vector<std::complex<double>>> {
 public:
  using DealiiVector = dealii::Vector<double>;
  using ComplexVector = std::vector<std::complex<double>>;
  ComplexVector Convert(const DealiiVector &input) const override {
    return ComplexVector{};
  };
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_
