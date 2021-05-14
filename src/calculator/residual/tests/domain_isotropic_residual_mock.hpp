#ifndef BART_SRC_CALCULATOR_RESIDUAL_TESTS_DOMAIN_ISOTROPIC_RESIDUAL_MOCK_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_TESTS_DOMAIN_ISOTROPIC_RESIDUAL_MOCK_HPP_

#include "calculator/residual/domain_isotropic_residual_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::residual {

class DomainIsotropicResidualMock : public DomainIsotropicResidualI {
 public:
  using Vector = dealii::Vector<double>;
  using FluxMoments = system::moments::SphericalHarmonicI;

  MOCK_METHOD(Vector, CalculateDomainResidual, (FluxMoments *, FluxMoments *), (override));
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_TESTS_DOMAIN_ISOTROPIC_RESIDUAL_MOCK_HPP_
