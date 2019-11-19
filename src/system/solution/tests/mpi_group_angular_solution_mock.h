#ifndef BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_
#define BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_

#include "system/solution/mpi_group_angular_solution_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system {

namespace solution {

class MPIGroupAngularSolutionMock : public MPIGroupAngularSolutionI {
 public:
  MOCK_CONST_METHOD0(total_angles, int());
  MOCK_CONST_METHOD0(solutions, SolutionMap&());
  MOCK_CONST_METHOD1(BracketOp, const MPIVector&(const AngleIndex));
  MOCK_METHOD1(BracketOp, MPIVector&(const AngleIndex));
  MOCK_METHOD1(GetSolution, MPIVector&(const AngleIndex));

  virtual const MPIVector& operator[](const AngleIndex angle) const override {
    return BracketOp(angle);};

  virtual MPIVector& operator[](const AngleIndex angle) override {
    return BracketOp(angle);};
};

} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_