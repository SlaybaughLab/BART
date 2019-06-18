#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "system/system_types.h"
#include "system/moments/spherical_harmonic_types.h"
#include "quadrature/angular/tests/angular_quadrature_set_mock.h"
#include "system/solution/tests/mpi_angular_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::Ref, ::testing::Return, ::testing::ReturnRef;

void SetVector(system::MPIVector& to_set, double value) {
  auto [first_row, last_row] = to_set.local_range();
  for (unsigned int i = first_row; i < last_row; ++i)
    to_set[i] = value;
  to_set.compress(dealii::VectorOperation::insert);
}

template <typename DimensionWrapper>
class QuadCalcSphericalHarmonicMomentsOnlyScalar : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  // Aliases
  using AngularQuadratureSetType = quadrature::angular::AngularQuadratureSetMock<dim>;
  using MomentCalculatorType = quadrature::calculators::SphericalHarmonicZerothMoment<dim>;

  // Pointer to tested object
  std::unique_ptr<MomentCalculatorType> test_calculator;

  // Supporting objects
  system::solution::MPIAngularMock mock_solution_;
  std::array<system::MPIVector, 3> mpi_vectors_;

  // Test object dependency
  std::shared_ptr<AngularQuadratureSetType> mock_angular_quad_;

  // Observing pointers
  AngularQuadratureSetType* angular_quad_obs_ptr_;

  // Test parameters
  const int n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int n_entries_per_proc = 10;

  void SetUp() override;
};

template <typename DimensionWrapper>
void QuadCalcSphericalHarmonicMomentsOnlyScalar<DimensionWrapper>::SetUp() {

  // Instantiate mock objects
  mock_angular_quad_ = std::make_shared<AngularQuadratureSetType>();

  // Instantiate object to be tested
  test_calculator = std::make_unique<MomentCalculatorType>(mock_angular_quad_);

  // Set up observation pointers
  angular_quad_obs_ptr_ = dynamic_cast<AngularQuadratureSetType*>(
      test_calculator->angular_quadrature_set_ptr());


  for (auto& mpi_vector : mpi_vectors_) {
    mpi_vector.reinit(MPI_COMM_WORLD,
                      n_processes * n_entries_per_proc,
                      n_entries_per_proc);
  };

  SetVector(mpi_vectors_[0], 1);
  SetVector(mpi_vectors_[1], 10);
  SetVector(mpi_vectors_[2], 100);
}

TYPED_TEST_CASE(QuadCalcSphericalHarmonicMomentsOnlyScalar,
                bart::testing::AllDimensions);

TYPED_TEST(QuadCalcSphericalHarmonicMomentsOnlyScalar, Constructor) {
  auto& mock_angular_quad_ = this->mock_angular_quad_;
  auto& angular_quad_obs_ptr_ = this->angular_quad_obs_ptr_;

  EXPECT_EQ(mock_angular_quad_.use_count(), 2);
  EXPECT_THAT(*mock_angular_quad_, Ref(*angular_quad_obs_ptr_));
}

TYPED_TEST(QuadCalcSphericalHarmonicMomentsOnlyScalar, CalculateBandAngleNumber) {
  auto& angular_quad_mock = *(this->angular_quad_obs_ptr_);
  auto& test_calculator = this->test_calculator;
  auto& mock_solution = this->mock_solution_;

  EXPECT_CALL(angular_quad_mock, total_quadrature_points())
      .WillOnce(Return(4));
  EXPECT_CALL(mock_solution, total_angles())
      .WillOnce(Return(3));
  EXPECT_ANY_THROW(test_calculator->CalculateMoment(&mock_solution, 0, 0, 0));
}

TYPED_TEST(QuadCalcSphericalHarmonicMomentsOnlyScalar, CalculateMomentsMPI) {
  auto& angular_quad_mock = *(this->angular_quad_obs_ptr_);
  auto& test_calculator = this->test_calculator;
  auto& mock_solution = this->mock_solution_;

  const int n_angles = 3;
  const int group = 0;

  EXPECT_CALL(angular_quad_mock, total_quadrature_points())
      .WillOnce(Return(n_angles));
  EXPECT_CALL(mock_solution, total_angles())
      .WillOnce(Return(n_angles));

  for (int angle = 0; angle < n_angles; ++angle) {
    system::Index index{group, angle};
    EXPECT_CALL(mock_solution, BracketOp(index))
        .WillOnce(ReturnRef(this->mpi_vectors_[angle]));
  }

  std::vector<quadrature::angular::Weight> weights{2.2, 3.3, 4.4};
  EXPECT_CALL(angular_quad_mock, quadrature_weights())
      .WillOnce(Return(weights));

  system::moments::MomentVector expected_result(
      this->n_entries_per_proc*this->n_processes);

  expected_result = 4.4*100 + 3.3*10 + 2.2;

  auto result = test_calculator->CalculateMoment(&mock_solution, group, 0, 0);
  EXPECT_EQ(result, expected_result);
}


} // namespace
