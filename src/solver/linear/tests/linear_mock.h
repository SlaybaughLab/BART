#ifndef BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_
#define BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_

#include "solver/linear/linear_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace solver {

class LinearMock : public LinearI {
 public:
  MOCK_METHOD(void, Solve, (dealii::PETScWrappers::MatrixBase *,
      dealii::PETScWrappers::VectorBase *,
      dealii::PETScWrappers::VectorBase *,
      dealii::PETScWrappers::PreconditionerBase *),
              (override));
};

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_
