#ifndef BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_
#define BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_

#include "formulation/stamper_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

template <int dim>
class StamperMock : public StamperI<dim> {
 public:
  MOCK_METHOD(void,
              StampMatrix,
              (system::MPISparseMatrix& to_stamp,
               std::function<void(formulation::FullMatrix&,
                                  const domain::CellPtr<dim>&)> stamp_function),
              (override));
  MOCK_METHOD(void,
              StampVector,
              (system::MPIVector& to_stamp,
              std::function<void(formulation::Vector&,
                                 const domain::CellPtr<dim>&)> stamp_function),
              (override));
  MOCK_METHOD(void,
              StampBoundaryMatrix,
              (system::MPISparseMatrix& to_stamp,
                  std::function<void(formulation::FullMatrix&,
                                     const domain::FaceIndex,
                                     const domain::CellPtr<dim>&)> stamp_function),
              (override));
  MOCK_METHOD(void,
              StampBoundaryVector,
              (system::MPIVector& to_stamp,
                  std::function<void(formulation::Vector&,
                                     const domain::FaceIndex,
                                     const domain::CellPtr<dim>&)> stamp_function),
              (override));
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_TESTS_STAMPER_MOCK_HPP_
