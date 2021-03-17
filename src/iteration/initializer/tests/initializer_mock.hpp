#ifndef BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_HPP_
#define BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_HPP_

#include "iteration/initializer/initializer_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart::iteration::initializer {

class InitializerMock : public InitializerI {
 public:
  MOCK_METHOD(void, Initialize, (system::System&), (override));
};

} // namespace bart::iteration::initializer

#endif //BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_HPP_
