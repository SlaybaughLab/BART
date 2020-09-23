#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_CALCULATOR_VECTOR_SUBTRACTOR_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_CALCULATOR_VECTOR_SUBTRACTOR_H_

#include <deal.II/lac/vector.h>

#include "instrumentation/converter/converter_i.h"
#include "utility/named_type.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace calculator {

using AbsoluteValue = utility::NamedType<bool, struct AbsoluteValueStruct>;

/*! \brief This converter returns the result of subtracting the provided vector
 *         from a stored vector.
 *
 *  This can be used to, for example, calculate the error in a solve. Options
 *  are provided to return the absolute value of the subtraction.
 */
class VectorSubtractor : public ConverterI<dealii::Vector<double>,
                                           dealii::Vector<double>> {
 public:
  using DealiiVector = dealii::Vector<double>;

  VectorSubtractor(DealiiVector minuend, AbsoluteValue calculate_absolute_value)
      : minuend_(minuend),
        calculate_absolute_value_(calculate_absolute_value) {}
  DealiiVector Convert(const DealiiVector &input) const override {
    return DealiiVector();
  }

  DealiiVector minuend() const { return minuend_; }
  bool calculate_absolute_value() const { return calculate_absolute_value_; }

 private:
  DealiiVector minuend_;
  bool calculate_absolute_value_{false};
};

} // namespace calculator

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_CALCULATOR_VECTOR_SUBTRACTOR_H_
