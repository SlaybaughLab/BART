#ifndef BART_SRC_CALCULATOR_RESIDUAL_TESTS_VECTOR_DIFFERENCE_MOCK_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_TESTS_VECTOR_DIFFERENCE_MOCK_HPP_

#include "calculator/residual/vector_difference_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::residual {

class VectorDifferenceMock : public VectorDifferenceI {
 public:
  MOCK_METHOD(Vector, CalculateResidual, (const Vector &, const Vector &), (const, override));
  MOCK_METHOD(Vector, CalculateResidual, (const Vector &, const Vector &, const double weight), (const, override));
};


} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_TESTS_VECTOR_DIFFERENCE_MOCK_HPP_
