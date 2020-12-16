#ifndef BART_SRC_ITERATION_SUBROUTINE_TESTS_SUBROUTINE_MOCK_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_TESTS_SUBROUTINE_MOCK_HPP_

#include "iteration/subroutine/subroutine_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::iteration::subroutine {

class SubroutineMock : public SubroutineI {
  MOCK_METHOD(void, Execute, (system::System&), (override));
};

} // namespace bart::iteration::subroutine


#endif //BART_SRC_ITERATION_SUBROUTINE_TESTS_SUBROUTINE_MOCK_HPP_
