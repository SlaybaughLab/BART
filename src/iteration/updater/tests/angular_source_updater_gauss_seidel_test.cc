#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>
#include <iteration/updater/angular_source_updater_gauss_seidel.h>
#include <deal.II/base/mpi.h>

#include "test_helpers/gmock_wrapper.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "formulation/tests/angular_stamper_mock.h"
#include "system/terms/tests/linear_term_mock.h"

namespace  {

using namespace bart;
using ::testing::A, ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class IterationUpdaterAngularSourceUpdaterGaussSeidelTest :
    public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using RightHandSideType = bart::system::terms::LinearTermMock;

  // Tested object
  std::unique_ptr<UpdaterType> test_updater_;

  // Dependencies
  std::shared_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;

  // Other test objects
  bart::system::System test_system_;
  std::shared_ptr<system::MPIVector> source_vector_ptr_;
  bart::system::MPIVector expected_vector_;
  RightHandSideType* right_hand_side_obs_ptr_;

  void SetUp() override;
  void SetUpTestObject();
  void SetUpSystem();
};

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUp() {
  SetUpTestObject();
  SetUpSystem();
}

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUpTestObject() {
  stamper_ptr_ = std::make_shared<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  test_updater_ = std::make_unique<UpdaterType>(stamper_ptr_, quadrature_set_ptr_);
}

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUpSystem() {
  auto mock_right_hand_side_ptr = std::make_unique<RightHandSideType>();
  right_hand_side_obs_ptr_ = mock_right_hand_side_ptr.get();

  source_vector_ptr_ = std::make_shared<system::MPIVector>();
  auto n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  source_vector_ptr_->reinit(MPI_COMM_WORLD, n_processes*5, 5);
  expected_vector_.reinit(*source_vector_ptr_);

  ON_CALL(*mock_right_hand_side_ptr, GetVariableTermPtr(A<system::Index>(), _))
      .WillByDefault(Return(source_vector_ptr_));

  test_system_.right_hand_side_ptr_ = std::move(mock_right_hand_side_ptr);
}



TYPED_TEST_SUITE(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest, Constructor) {
  constexpr int dim = this->dim;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  EXPECT_NO_THROW({ UpdaterType test_updater(stamper_ptr, quadrature_set_ptr); });
  EXPECT_ANY_THROW({ UpdaterType test_updater(nullptr, quadrature_set_ptr); });
  EXPECT_ANY_THROW({ UpdaterType test_updater(stamper_ptr, nullptr); });
}

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest, Getters) {
  auto quadrature_set_ptr = this->test_updater_->quadrature_set_ptr();
  auto stamper_ptr = this->test_updater_->stamper_ptr();

  EXPECT_EQ(quadrature_set_ptr, this->quadrature_set_ptr_.get());
  EXPECT_EQ(stamper_ptr, this->stamper_ptr_.get());
}



} // namespace
