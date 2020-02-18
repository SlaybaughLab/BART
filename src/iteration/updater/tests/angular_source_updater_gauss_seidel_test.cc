#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>
#include <iteration/updater/angular_source_updater_gauss_seidel.h>
#include <deal.II/base/mpi.h>

#include "test_helpers/gmock_wrapper.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "formulation/tests/angular_stamper_mock.h"
#include "system/terms/tests/linear_term_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.h"

namespace  {

using namespace bart;
using ::testing::A, ::testing::An, ::testing::Return, ::testing::_;
using ::testing::NiceMock, ::testing::DoDefault, ::testing::WithArg;
using ::testing::Invoke, ::testing::Ref, ::testing::ReturnRef;

template <typename DimensionWrapper>
class IterationUpdaterAngularSourceUpdaterGaussSeidelTest :
    public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;
  using QuadraturePointType = quadrature::QuadraturePointI<dim>;
  using RightHandSideType = NiceMock<bart::system::terms::LinearTermMock>;
  using MomentsType = bart::system::moments::SphericalHarmonicMock;

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
  MomentsType* current_moments_obs_ptr_;
  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;
  bart::system::moments::MomentsMap current_iteration_moments_,
      previous_iteration_moments_;

  // Test parameters
  const int total_groups = test_helpers::RandomDouble(2, 5);
  const int total_angles = test_helpers::RandomDouble(2, 5);
  const int l_max = 2;
  const int group = test_helpers::RandomDouble(0, total_groups + 1); // group number
  const int angle_number = test_helpers::RandomDouble(0, total_angles + 1);
  const bart::system::Index solution_index{group, angle_number};

  void SetUp() override;
  void SetUpTestObject();
  void SetUpSystem();
};

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUp() {
  SetUpTestObject();
  SetUpSystem();
  quadrature_point_ptr_ = std::make_shared<quadrature::QuadraturePointMock<dim>>();
}

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUpTestObject() {
  stamper_ptr_ = std::make_shared<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  test_updater_ = std::make_unique<UpdaterType>(stamper_ptr_, quadrature_set_ptr_);
}

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUpSystem() {
  // RHS MPI Vectors
  auto mock_right_hand_side_ptr = std::make_unique<RightHandSideType>();
  right_hand_side_obs_ptr_ = mock_right_hand_side_ptr.get();

  source_vector_ptr_ = std::make_shared<system::MPIVector>();
  auto n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  source_vector_ptr_->reinit(MPI_COMM_WORLD, n_processes*5, 5);
  expected_vector_.reinit(*source_vector_ptr_);

  ON_CALL(*mock_right_hand_side_ptr, GetVariableTermPtr(A<system::Index>(), _))
      .WillByDefault(Return(source_vector_ptr_));

  test_system_.right_hand_side_ptr_ = std::move(mock_right_hand_side_ptr);

  // Moments
  auto current_moments_ptr_ = std::make_unique<system::moments::SphericalHarmonicMock>();
  current_moments_obs_ptr_ = current_moments_ptr_.get();

  for (system::GroupNumber group = 0; group < total_groups; ++group) {
    for (system::moments::HarmonicL l = 0; l < l_max; ++l) {
      for (system::moments::HarmonicM m = -l_max; m <= l_max; ++m) {
        system::moments::MomentVector current_moment, previous_moment;
        current_iteration_moments_[{group, l, m}] = current_moment;
        previous_iteration_moments_[{group, l, m}] = previous_moment;
      }
    }
  }

  test_system_.current_moments = std::move(current_moments_ptr_);
}

// ===== HELPER FUNCTIONS ======================================================

// Fills an MPI vector with value
void StampMPIVector(bart::system::MPIVector &to_fill, double value = 2) {
  auto [local_begin, local_end] = to_fill.local_range();
  for (unsigned int i = local_begin; i < local_end; ++i)
    to_fill(i) += value;
  to_fill.compress(dealii::VectorOperation::add);
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

/*
 * ======== UpdateScatteringSource Tests =======================================
 */

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
    UpdateScatteringSourceBadRHS) {
  using term = bart::system::terms::VariableLinearTerms;
  EXPECT_CALL(*this->right_hand_side_obs_ptr_,
              GetVariableTermPtr(this->solution_index, term::kScatteringSource))
      .WillOnce(Return(nullptr));
  EXPECT_ANY_THROW({
    this->test_updater_->UpdateScatteringSource(this->test_system_,
                                                this->group,
                                                this->angle_number);
  });
}

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
    UpdateScatteringSource) {
  using term = bart::system::terms::VariableLinearTerms;
  StampMPIVector(*this->source_vector_ptr_, 3); // Fill with a random value, should be zero'd
  double expected_value = test_helpers::RandomDouble(1, 20);
  StampMPIVector(this->expected_vector_, expected_value);

  auto stamper_function = [&](bart::system::MPIVector &to_stamp) {
    StampMPIVector(to_stamp, expected_value);
  };

  EXPECT_CALL(*this->right_hand_side_obs_ptr_,
              GetVariableTermPtr(this->solution_index, term::kScatteringSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_,
              GetQuadraturePoint(quadrature::QuadraturePointIndex(this->angle_number)))
      .WillOnce(Return(this->quadrature_point_ptr_));
  EXPECT_CALL(*this->stamper_ptr_,
              StampScatteringSourceTerm(Ref(*this->source_vector_ptr_),
                  this->quadrature_point_ptr_,
                  system::EnergyGroup(this->group),
                  Ref(this->current_iteration_moments_[{this->group, 0, 0}]),
                  Ref(this->current_iteration_moments_)))
      .WillOnce(WithArg<0>(Invoke(stamper_function)));
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(ReturnRef(this->current_iteration_moments_));

  this->test_updater_->UpdateScatteringSource(this->test_system_,
                                              this->group,
                                              this->angle_number);
  EXPECT_TRUE(bart::testing::CompareMPIVectors(this->expected_vector_,
                                               *this->source_vector_ptr_));
}

/*
 * ======== UpdateFissionSource Tests ==========================================
 */

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
           UpdateFissionSourceBadRHS) {
  using term = bart::system::terms::VariableLinearTerms;
  this->test_system_.k_effective = 1.0;
  EXPECT_CALL(*this->right_hand_side_obs_ptr_,
              GetVariableTermPtr(this->solution_index, term::kFissionSource))
      .WillOnce(Return(nullptr));
  EXPECT_ANY_THROW({
    this->test_updater_->UpdateFissionSource(this->test_system_,
                                             this->group,
                                             this->angle_number);
                   });
}

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
           UpdateFissionSource) {
  using term = bart::system::terms::VariableLinearTerms;
  const double k_eff = 1.1023;
  this->test_system_.k_effective = k_eff;
  StampMPIVector(*this->source_vector_ptr_, 3); // Fill with a value, should be zero'd
  double expected_value = test_helpers::RandomDouble(1, 20);
  StampMPIVector(this->expected_vector_, expected_value);

  auto stamper_function = [&](bart::system::MPIVector &to_stamp) {
    StampMPIVector(to_stamp, expected_value);
  };

  EXPECT_CALL(*this->right_hand_side_obs_ptr_,
              GetVariableTermPtr(this->solution_index, term::kFissionSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_,
              GetQuadraturePoint(quadrature::QuadraturePointIndex(this->angle_number)))
      .WillOnce(Return(this->quadrature_point_ptr_));
  EXPECT_CALL(*this->stamper_ptr_,
              StampFissionSourceTerm(Ref(*this->source_vector_ptr_),
                                     this->quadrature_point_ptr_,
                                     system::EnergyGroup(this->group),
                                     k_eff,
                                     Ref(this->current_iteration_moments_[{this->group, 0, 0}]),
                                     Ref(this->current_iteration_moments_)))
      .WillOnce(WithArg<0>(Invoke(stamper_function)));
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(ReturnRef(this->current_iteration_moments_));

  this->test_updater_->UpdateFissionSource(this->test_system_,
                                           this->group,
                                           this->angle_number);
  EXPECT_TRUE(bart::testing::CompareMPIVectors(this->expected_vector_,
                                               *this->source_vector_ptr_));
}

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
           UpdateFissionSourceBadKeff) {
  using term = bart::system::terms::VariableLinearTerms;

  EXPECT_CALL(*this->right_hand_side_obs_ptr_,
              GetVariableTermPtr(this->solution_index, term::kFissionSource))
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_,
              GetQuadraturePoint(quadrature::QuadraturePointIndex(this->angle_number)))
      .WillRepeatedly(Return(this->quadrature_point_ptr_));
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillRepeatedly(ReturnRef(this->current_iteration_moments_));

  std::vector<std::optional<double>> bad_k_eff{-1.0, 0.0, std::nullopt};

  for (const auto k_eff : bad_k_eff) {
    this->test_system_.k_effective = k_eff;
    EXPECT_ANY_THROW({
    this->test_updater_->UpdateFissionSource(this->test_system_,
                                             this->group,
                                             this->angle_number);
                     });
  }
}

} // namespace
