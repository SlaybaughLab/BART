#ifndef BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_HPP_

#include "calculator/residual/vector_difference_i.hpp"

namespace bart::calculator::residual {

/*! \brief Default implementation for VectorDifference.
 *
 * This acts almost singly as a wrapper of dealii functions.
 *
 */
class VectorDifference : public VectorDifferenceI {
 public:
  virtual ~VectorDifference() = default;
  [[nodiscard]] auto CalculateResidual(const Vector& minuend, const Vector& subtrahend) const -> Vector override {
    return CalculateResidual(minuend, subtrahend, 1);
  }

  [[nodiscard]] auto CalculateResidual(const Vector& minuend, const Vector& subtrahend,
                                       const double weight) const -> Vector override {
    Vector return_value = minuend;
    return_value.add(-1, subtrahend);
    return_value *= weight;
    return return_value;
  }
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_HPP_
