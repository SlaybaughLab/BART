#ifndef BART_SRC_EIGENVALUE_K_EFFECTIVE_TESTS_K_EFFECTIVE_UPDATER_MOCK_H_
#define BART_SRC_EIGENVALUE_K_EFFECTIVE_TESTS_K_EFFECTIVE_UPDATER_MOCK_H_

#include "eigenvalue/k_effective/k_effective_updater_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace eigenvalue {

namespace k_effective {

class K_EffectiveUpdaterMock : public K_EffectiveUpdaterI {
 public:
  MOCK_CONST_METHOD0(k_effective, std::optional<double>());
  MOCK_METHOD1(CalculateK_Effective, double(system::System& system));
};

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart

#endif //BART_SRC_EIGENVALUE_K_EFFECTIVE_TESTS_K_EFFECTIVE_UPDATER_MOCK_H_
