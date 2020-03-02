#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_

#include <functional>
#include <memory>

#include "domain/domain_types.h"
#include "formulation/formulation_types.h"
#include "formulation/tests/stamper_mock.h"
#include "system/terms/tests/bilinear_term_mock.h"
#include "system/terms/tests/linear_term_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/dealii_test_domain.h"

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

  // Test parameters
  const int total_groups = bart::test_helpers::RandomDouble(2, 5);

  // Call parameters (objects and parameters used in update calls)
  system::System test_system_;
  const int group_number = bart::test_helpers::RandomDouble(0, total_groups);
  const int angle_index = bart::test_helpers::RandomDouble(0, 10);
  const system::Index index{group_number, angle_index};

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
  static void EvaluateFunction(std::function<void(formulation::FullMatrix&,
                                                  const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateVectorFunction(std::function<void(formulation::Vector&,
                                                        const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateBoundaryFunction(std::function<void(formulation::FullMatrix&,
                                                          const domain::FaceIndex,
                                                          const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateVectorBoundaryFunction(std::function<void(formulation::Vector&,
                                                                const domain::FaceIndex,
                                                                const domain::CellPtr<dim>&)> stamp_function);
 private:
  void SetUpSystem();
};

template <int dim>
void UpdaterTests<dim>::SetUp() {
  this->SetUpDealii();
  this->SetUpSystem();
}

template <int dim>
std::unique_ptr<StamperMock<dim>> UpdaterTests<dim>::MakeStamper() {
  auto mock_stamper_ptr = std::make_unique<StamperMock<dim>>();

  ON_CALL(*mock_stamper_ptr, StampMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateFunction)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateBoundaryFunction)));
  ON_CALL(*mock_stamper_ptr, StampVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateVectorFunction)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateVectorBoundaryFunction)));

  return std::move(mock_stamper_ptr);
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
  ON_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(A<system::Index>()))
      .WillByDefault(Return(this->matrix_to_stamp));
  ON_CALL(*current_moments_obs_ptr_, moments())
      .WillByDefault(ReturnRef(current_iteration_moments_));
}

template <int dim>
void UpdaterTests<dim>::EvaluateFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  stamp_function(to_stamp, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorFunction(
    std::function<void(formulation::Vector&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  domain::CellPtr<dim> cell_ptr;
  stamp_function(to_stamp, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateBoundaryFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  domain::FaceIndex index(0);
  stamp_function(to_stamp, index, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorBoundaryFunction(
    std::function<void(formulation::Vector&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  domain::CellPtr<dim> cell_ptr;
  domain::FaceIndex index(0);
  stamp_function(to_stamp, index, cell_ptr);
}

} // namespace test_helpers

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_
