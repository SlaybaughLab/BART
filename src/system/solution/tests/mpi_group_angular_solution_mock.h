#ifndef BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_
#define BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_

#include "system/solution/mpi_group_angular_solution_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system {

namespace solution {

class MPIGroupAngularSolutionMock : public MPIGroupAngularSolutionI {
 public:
  MOCK_METHOD(int, total_angles, (), (override, const));
  MOCK_METHOD(const SolutionMap&, solutions, (), (override, const));
  MOCK_METHOD(SolutionMap&, solutions, (), (override));
  MOCK_METHOD(const MPIVector&, BracketOp, (const AngleIndex), (const));
  MOCK_METHOD(MPIVector&, BracketOp, (const AngleIndex));
  MOCK_METHOD(MPIVector&, GetSolution, (const AngleIndex), (override));

  virtual const MPIVector& operator[](const AngleIndex angle) const override {
    return BracketOp(angle);};

  virtual MPIVector& operator[](const AngleIndex angle) override {
    return BracketOp(angle);};
};

} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_GROUP_ANGULAR_SOLUTION_MOCK_H_