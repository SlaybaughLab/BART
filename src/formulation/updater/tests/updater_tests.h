#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_

#include <functional>
#include <memory>

#include "domain/domain_types.hpp"
#include "formulation/formulation_types.hpp"
#include "formulation/tests/stamper_mock.hpp"
#include "system/terms/tests/bilinear_term_mock.hpp"
#include "system/terms/tests/linear_term_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/system.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/dealii_test_domain.h"
#include "problem/parameter_types.hpp"

namespace bart {

namespace formulation {

namespace updater {

namespace test_helpers {

using ::testing::WithArg, ::testing::Invoke, ::testing::_, ::testing::ReturnRef,
::testing::Return, ::testing::A;

template <int dim>
class UpdaterTests : public ::testing::Test,
                     public bart::testing::DealiiTestDomain<dim> {
 public:
  using StamperType = formulation::StamperMock<dim>;
  using MomentsType = system::moments::SphericalHarmonicMock;
  using SystemType = system::System;

  // Test parameters
  const int total_groups = bart::test_helpers::RandomDouble(2, 5);
  const int total_angles = total_groups + 1;

  // Call parameters (objects and parameters used in update calls)
  SystemType test_system_;
  const int group_number = bart::test_helpers::RandomDouble(0, total_groups);
  const int angle_index = bart::test_helpers::RandomDouble(0, 10);
  const int reflected_angle_index = bart::test_helpers::RandomDouble(11, 20);
  const system::Index index{group_number, angle_index};
  const system::Index reflected_index{group_number, reflected_angle_index};

  // Mock system object observation pointers
  bart::system::terms::BilinearTermMock* mock_lhs_obs_ptr_;
  bart::system::terms::LinearTermMock* mock_rhs_obs_ptr_;
  MomentsType* current_moments_obs_ptr_;

  // Real system objects
  bart::system::moments::MomentsMap current_iteration_moments_;

  // Vectors to "stamp", these will be filled with random values that should be
  // properly zero'd
  std::shared_ptr<bart::system::MPISparseMatrix> matrix_to_stamp;
  std::shared_ptr<bart::system::MPIVector> vector_to_stamp;
  // Vectors to compare the stamped vectors to, these will be set to 0
  system::MPISparseMatrix& expected_result = this->matrix_2;
  system::MPIVector& expected_vector_result = this->vector_2;

  // Calls To be used by derived tests
  std::unique_ptr<StamperMock<dim>> MakeStamper();
  void SetUp() override;

  // Functions invoked by the Stamper to actually evaluate the passed function
  void EvaluateMatrixFunctionOnDomain(std::function<void(formulation::FullMatrix&,
                                                         const domain::CellPtr<dim>&)> stamp_function);
  void EvaluateVectorFunctionOnDomain(std::function<void(formulation::Vector&,
                                                         const domain::CellPtr<dim>&)> stamp_function);
  void EvaluateMatrixFunctionOnBoundary(std::function<void(formulation::FullMatrix&,
                                                          const domain::FaceIndex,
                                                          const domain::CellPtr<dim>&)> stamp_function);
  void EvaluateVectorFunctionOnBoundary(std::function<void(formulation::Vector&,
                                                           const domain::FaceIndex,
                                                           const domain::CellPtr<dim>&)> stamp_function);
 private:
  void SetUpSystem();
  void SetUpBoundaries();
};

template <int dim>
void UpdaterTests<dim>::SetUp() {
  this->SetUpDealii();
  this->SetUpSystem();
  this->SetUpBoundaries();
}

template <int dim>
std::unique_ptr<StamperMock<dim>> UpdaterTests<dim>::MakeStamper() {
  auto mock_stamper_ptr = std::make_unique<StamperMock<dim>>();

  ON_CALL(*mock_stamper_ptr, StampMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this, &formulation::updater::test_helpers::UpdaterTests<dim>::EvaluateMatrixFunctionOnDomain)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this, &formulation::updater::test_helpers::UpdaterTests<dim>::EvaluateMatrixFunctionOnBoundary)));
  ON_CALL(*mock_stamper_ptr, StampVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this, &formulation::updater::test_helpers::UpdaterTests<dim>::EvaluateVectorFunctionOnDomain)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this, &formulation::updater::test_helpers::UpdaterTests<dim>::EvaluateVectorFunctionOnBoundary)));

  return mock_stamper_ptr;
}

template <int dim>
void UpdaterTests<dim>::SetUpSystem() {
  auto mock_lhs_ptr = std::make_unique<system::terms::BilinearTermMock>();
  mock_lhs_obs_ptr_ = mock_lhs_ptr.get();
  auto mock_rhs_ptr = std::make_unique<system::terms::LinearTermMock>();
  mock_rhs_obs_ptr_ = mock_rhs_ptr.get();

  auto current_moments_ptr = std::make_unique<system::moments::SphericalHarmonicMock>();
  current_moments_obs_ptr_ = current_moments_ptr.get();

  for (int group = 0; group < this->total_groups; ++group) {
    system::moments::MomentVector new_vector;
    current_iteration_moments_.insert_or_assign({group, 0, 0}, new_vector);
  }

  matrix_to_stamp = std::make_shared<system::MPISparseMatrix>();
  matrix_to_stamp->reinit(this->matrix_1);
  vector_to_stamp = std::make_shared<system::MPIVector>();
  vector_to_stamp->reinit(this->vector_1);
  // Set matrix, vector to stamp to a random value
  double random_value = bart::test_helpers::RandomDouble(1, 10);
  this->StampMatrix(*this->matrix_to_stamp, random_value);
  auto local_elements = this->vector_to_stamp->locally_owned_elements();
  for (auto entry : local_elements)
    this->vector_to_stamp->operator()(entry) = random_value;
  this->vector_to_stamp->compress(dealii::VectorOperation::insert);
  // Set expected vectors and matrices to zero
  expected_vector_result = 0;
  expected_result = 0;

  test_system_.left_hand_side_ptr_ = std::move(mock_lhs_ptr);
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr);
  test_system_.current_moments = std::move(current_moments_ptr);

  ON_CALL(*this->mock_rhs_obs_ptr_, GetVariableTermPtr(A<system::Index>(), _))
      .WillByDefault(Return(this->vector_to_stamp));
  ON_CALL(*this->mock_rhs_obs_ptr_, GetFixedTermPtr(A<system::Index>()))
      .WillByDefault(Return(this->vector_to_stamp));
  ON_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(A<system::Index>()))
      .WillByDefault(Return(this->matrix_to_stamp));
  ON_CALL(*current_moments_obs_ptr_, moments())
      .WillByDefault(ReturnRef(current_iteration_moments_));
}

template <int dim>
void UpdaterTests<dim>::SetUpBoundaries() {
  using Boundary = bart::problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  double zero_tol = 1.0e-14;

  for (auto &cell : this->cells_) {
    for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
      auto face = cell->face(face_id);
      dealii::Point<dim> face_center = face->center();

      switch (dim) {
        case 3: {
          if (std::fabs(face_center[2]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kZMin));
            break;
          } else if (std::fabs(face_center[2] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kZMax));
            break;
          }
          [[fallthrough]];
        }
        case 2: {
          if (std::fabs(face_center[1]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kYMin));
            break;
          } else if (std::fabs(face_center[1] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kYMax));
            break;
          }
          [[fallthrough]];
        }
        case 1: {
          if (std::fabs(face_center[0]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kXMin));
            break;
          } else if (std::fabs(face_center[0] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kXMax));
            break;
          }
        }
      }
    }
  }
}

template <int dim>
void UpdaterTests<dim>::EvaluateMatrixFunctionOnDomain(
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  for (auto& cell_ptr : this->cells_) {
    stamp_function(to_stamp, cell_ptr);
  }
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorFunctionOnDomain(
    std::function<void(formulation::Vector&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  for (auto& cell_ptr : this->cells_) {
    stamp_function(to_stamp, cell_ptr);
  }
}

template <int dim>
void UpdaterTests<dim>::EvaluateMatrixFunctionOnBoundary(
    std::function<void(formulation::FullMatrix&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  for (auto& cell_ptr : this->cells_) {
    if (cell_ptr->at_boundary()) {
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell_ptr->face(face)->at_boundary()) {
          stamp_function(to_stamp, domain::FaceIndex(face), cell_ptr);
        }
      }
    }
  }
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorFunctionOnBoundary(
    std::function<void(formulation::Vector&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  for (auto& cell_ptr : this->cells_) {
    if (cell_ptr->at_boundary()) {
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell_ptr->face(face)->at_boundary()) {
          stamp_function(to_stamp, domain::FaceIndex(face), cell_ptr);
        }
      }
    }
  }
}

} // namespace test_helpers

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_
