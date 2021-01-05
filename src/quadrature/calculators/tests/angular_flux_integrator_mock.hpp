#ifndef BART_SRC_QUADRATURE_CALCULATORS_TESTS_ANGULAR_FLUX_INTEGRATOR_MOCK_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_TESTS_ANGULAR_FLUX_INTEGRATOR_MOCK_HPP_

#include "quadrature/calculators/angular_flux_integrator_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::quadrature::calculators {

class AngularFluxIntegratorMock : public AngularFluxIntegratorI {
 public:
  MOCK_METHOD(Vector, Integrate, (const VectorMap&), (const, override));
  MOCK_METHOD(Vector, NetCurrent, (const VectorMap&, const DegreeOfFreedom), (const, override));
};

} // namespace bart::quadrature::calculators



#endif //BART_SRC_QUADRATURE_CALCULATORS_TESTS_ANGULAR_FLUX_INTEGRATOR_MOCK_HPP_
