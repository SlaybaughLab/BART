#ifndef BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_
#define BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_

#include "formulation/stamper_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::formulation {

template <int dim>
class StamperMock : public StamperI<dim> {
 public:
  using typename StamperI<dim>::CellMatrixStampFunction;
  using typename StamperI<dim>::CellVectorStampFunction;
  using typename StamperI<dim>::FaceMatrixStampFunction;
  using typename StamperI<dim>::FaceVectorStampFunction;

  MOCK_METHOD(void, StampMatrix, (system::MPISparseMatrix&, CellMatrixStampFunction), (override));
  MOCK_METHOD(void, StampVector, (system::MPIVector& to_stamp, CellVectorStampFunction), (override));
  MOCK_METHOD(void, StampBoundaryMatrix, (system::MPISparseMatrix& to_stamp, FaceMatrixStampFunction), (override));
  MOCK_METHOD(void, StampBoundaryVector, (system::MPIVector& to_stamp, FaceVectorStampFunction), (override));
};

} // namespace bart::formulation

#endif //BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_
