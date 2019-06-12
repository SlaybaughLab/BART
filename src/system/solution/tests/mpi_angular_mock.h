#ifndef BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_ANGULAR_MOCK_H_
#define BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_ANGULAR_MOCK_H_

#include "system/solution/mpi_angular_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system {

namespace solution {

class MPIAngularMock : public MPIAngularI {
 public:
  MOCK_CONST_METHOD0(total_groups, int());
  MOCK_CONST_METHOD0(total_angles, int());
  MOCK_CONST_METHOD0(solutions, SolutionMap&());
  MOCK_CONST_METHOD1(BracketOp, const MPIVector&(const Index));
  MOCK_METHOD1(BracketOp, MPIVector&(const Index));

  virtual const MPIVector& operator[](const Index index) const override {
    return BracketOp(index);};

  virtual MPIVector& operator[](const Index index) override {
    return BracketOp(index);};
};

} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_TESTS_MPI_ANGULAR_MOCK_H_