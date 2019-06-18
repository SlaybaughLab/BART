#ifndef BART_SRC_ITERATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_
#define BART_SRC_ITERATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_

#include "iteration/updater/fixed_updater_i.h"

#include "system/system_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace data {

struct System;

} // namespace data

namespace iteration {

namespace updater {

class FixedUpdaterMock : public FixedUpdaterI {
 public:
  MOCK_METHOD3(UpdateFixedTerms, void(
      data::System& system,
      data::system::GroupNumber group,
      data::system::AngleIndex angle));
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_TESTS_FIXED_UPDATER_MOCK_H_