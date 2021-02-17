#ifndef BART_SRC_CALCULATOR_RESIDUAL_TESTS_ISOTROPIC_RESIDUAL_MOCK_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_TESTS_ISOTROPIC_RESIDUAL_MOCK_HPP_

#include "calculator/residual/isotropic_residual_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::residual {

class IsotropicResidualMock : public IsotropicResidualI {
 public:
  MOCK_METHOD(Vector, CalculateIsotropicResidual, (Moments *, Moments *, const int group, const FullMatrix &sigma_s),
              (const, override));
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_TESTS_ISOTROPIC_RESIDUAL_MOCK_HPP_
