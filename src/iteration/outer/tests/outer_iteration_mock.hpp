#ifndef BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_HPP_
#define BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_HPP_

#include "iteration/outer/outer_iteration_i.hpp"

#include "system/system.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterIterationMock : public OuterIterationI {
 public:
  using Subroutine = OuterIterationI::Subroutine;
  MOCK_METHOD(void, IterateToConvergence, (system::System &), (override));
  MOCK_METHOD(OuterIterationMock&, AddPostIterationSubroutine, (std::unique_ptr<Subroutine>), (override));
}; // namespace outer

} // namespace iteration


} // namespace system

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_TESTS_OUTER_ITERATION_MOCK_HPP_
