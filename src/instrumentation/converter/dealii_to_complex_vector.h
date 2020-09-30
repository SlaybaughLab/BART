#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

#include "instrumentation/converter/factory.h"
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
    const int vector_size = input.size();
    ComplexVector return_vector(vector_size);
    for (int i = 0; i < vector_size; ++i) {
      return_vector.at(i).real(input[i]);
      return_vector.at(i).imag(0);
    }
    return return_vector;
  };
 private:
  static inline bool is_registered_ =
      ConverterIFactory<DealiiVector, ComplexVector>::get()
      .RegisterConstructor(ConverterName::kDealiiToComplexVector,
          []() {
            std::unique_ptr<ConverterI<DealiiVector, ComplexVector>> return_ptr =
                std::make_unique<DealiiToComplexVector>();
            return return_ptr;
      });
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_DEALII_TO_COMPLEX_VECTOR_H_
