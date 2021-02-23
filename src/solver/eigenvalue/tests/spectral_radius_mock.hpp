#ifndef BART_SRC_SOLVER_EIGENVALUE_TESTS_SPECTRAL_RADIUS_MOCK_HPP_
#define BART_SRC_SOLVER_EIGENVALUE_TESTS_SPECTRAL_RADIUS_MOCK_HPP_

#include "solver/eigenvalue/spectral_radius_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::solver::eigenvalue {

class SpectralRadiusMock : public SpectralRadiusI {
 public:
  using SpectralRadiusI::MatrixBase;
  MOCK_METHOD((std::pair<double, std::vector<double>>), SpectralRadius, (const MatrixBase *base), (override));
};


} // namespace bart::solver::eigenvalue

#endif //BART_SRC_SOLVER_EIGENVALUE_TESTS_SPECTRAL_RADIUS_MOCK_HPP_
