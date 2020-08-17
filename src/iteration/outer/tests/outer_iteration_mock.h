#ifndef BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_H_
#define BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_H_

#include "iteration/outer/outer_iteration_i.h"

#include "system/system.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterIterationMock : public OuterIterationI {
 public:
  MOCK_METHOD(void, IterateToConvergence, (system::System &), (override));
}; // namespace outer

} // namespace iteration


} // namespace system

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_H_
