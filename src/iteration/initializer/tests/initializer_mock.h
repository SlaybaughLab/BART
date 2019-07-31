#ifndef BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_H_
#define BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_H_

#include "iteration/initializer/initializer_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace iteration {

namespace initializer {

class InitializerMock : public InitializerI {
 public:
  MOCK_METHOD1(Initialize, void(system::System& sys));
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_INITIALIZER_TESTS_INITIALIZER_MOCK_H_
