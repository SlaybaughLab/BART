#ifndef BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_
#define BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_

#include "solver/linear_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace solver {

class LinearMock : public LinearI {
 public:
  MOCK_METHOD4(Solve, void(
      dealii::PETScWrappers::MatrixBase *A,
      dealii::PETScWrappers::VectorBase *x,
      dealii::PETScWrappers::VectorBase *b,
      dealii::PETScWrappers::PreconditionerBase *preconditioner));
};

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_TESTS_LINEAR_MOCK_H_
