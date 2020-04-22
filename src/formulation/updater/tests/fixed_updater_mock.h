#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_

#include "formulation/updater/fixed_updater_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

namespace updater {

class FixedUpdaterMock : public FixedUpdaterI {
 public:
  MOCK_METHOD(void, UpdateFixedTerms,
      (system::System&, system::EnergyGroup, quadrature::QuadraturePointIndex),
      (override));
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_
