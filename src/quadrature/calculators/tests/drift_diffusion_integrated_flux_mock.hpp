#ifndef BART_SRC_QUADRATURE_CALCULATORS_TESTS_DRIFT_DIFFUSION_INTEGRATED_FLUX_MOCK_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_TESTS_DRIFT_DIFFUSION_INTEGRATED_FLUX_MOCK_HPP_

#include "quadrature/calculators/drift_diffusion_integrated_flux_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::quadrature::calculators {

class DriftDiffusionIntegratedFluxMock : public DriftDiffusionIntegratedFluxI {
 public:
  MOCK_METHOD(Vector, Integrate, (const VectorMap&), (const, override));
};

} // namespace bart::quadrature::calculators



#endif //BART_SRC_QUADRATURE_CALCULATORS_TESTS_DRIFT_DIFFUSION_INTEGRATED_FLUX_MOCK_HPP_
