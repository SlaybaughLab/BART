#ifndef BART_SRC_ITERATION_UPDATER_TESTS_SOURCE_UPDATER_MOCK_H_
#define BART_SRC_ITERATION_UPDATER_TESTS_SOURCE_UPDATER_MOCK_H_

#include "iteration/updater/source_updater_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace iteration {

namespace updater {

class SourceUpdaterMock : public SourceUpdaterI {
 public:
  MOCK_METHOD3(UpdateScatteringSource, void(system::System& system,
      system::GroupNumber group,
      system::AngleIndex angle));
  MOCK_METHOD3(UpdateFissionSource, void(system::System& system,
      system::GroupNumber group,
      system::AngleIndex angle));
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_UPDATER_TESTS_SOURCE_UPDATER_MOCK_H_
