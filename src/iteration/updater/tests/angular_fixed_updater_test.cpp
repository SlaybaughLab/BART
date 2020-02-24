#include <numeric>

#include "iteration/updater/angular_fixed_updater.h"
#include "formulation/tests/angular_stamper_mock.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "system/system.h"
#include "system/terms/tests/bilinear_term_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return, ::testing::Ref, ::testing::WithArg, ::testing::Invoke;

/*
 * ===== BASIC TESTS ===========================================================
 * Tests operation of the constructor, provides mock stamper and observing ptr
 * to stamper.
 *
 */

template <typename DimensionWrapper>
class IterationUpdaterAngularFixedUpdaterTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using FixedUpdaterType = iteration::updater::AngularFixedUpdater<formulation::AngularStamperI<dim>>;

  std::shared_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  std::unique_ptr<FixedUpdaterType> test_updater_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationUpdaterAngularFixedUpdaterTest<DimensionWrapper>::SetUp() {
  stamper_ptr_ = std::make_shared<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  test_updater_ptr_ = std::make_unique<FixedUpdaterType>(
      stamper_ptr_, quadrature_set_ptr_);
}

TYPED_TEST_SUITE(IterationUpdaterAngularFixedUpdaterTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularFixedUpdaterTest, Constructor) {
  constexpr int dim = this->dim;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using FixedUpdaterType =
      iteration::updater::AngularFixedUpdater<formulation::AngularStamperI<dim>>;

  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();

  EXPECT_NO_THROW({
    FixedUpdaterType test_updater(stamper_ptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, nullptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, quadrature_set_ptr);
  });
}

TYPED_TEST(IterationUpdaterAngularFixedUpdaterTest, Getters) {
  auto stamper_ptr = this->test_updater_ptr_->stamper_ptr();
  auto quadrature_set_ptr = this->test_updater_ptr_->quadrature_set_ptr();

  ASSERT_NE(nullptr, stamper_ptr);
  EXPECT_EQ(stamper_ptr, this->stamper_ptr_.get());

  ASSERT_NE(nullptr, quadrature_set_ptr);
  EXPECT_EQ(quadrature_set_ptr, this->quadrature_set_ptr_.get());
}

/* ===== DOMAIN TESTS ==========================================================
 * Tests updating of system matrices. DealiiTestDomain is used to provide
 * matrices with sparsity patterns that can be stamped properly. A test system
 * struct will be instantiated (populated with mock classes) and passed to
 * the updater to verify it mediates actions between the system and the
 * stamper correctly.
 */

template <typename DimensionWrapper>
class IterationUpdaterAngularFixedUpdaterDomainTest :
    public IterationUpdaterAngularFixedUpdaterTest<DimensionWrapper>,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  IterationUpdaterAngularFixedUpdaterDomainTest()
      : group_number_(test_helpers::RandomDouble(0, 10)),
        angle_index_(test_helpers::RandomDouble(0, 10)),
        index_({group_number_, angle_index_}) {};

  bart::system::System test_system_;
  std::shared_ptr<system::MPISparseMatrix> matrix_ptr_;
  bart::system::terms::BilinearTermMock* mock_lhs_obs_ptr_;

  const system::GroupNumber group_number_;
  const system::AngleIndex angle_index_;
  const system::Index index_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationUpdaterAngularFixedUpdaterDomainTest<DimensionWrapper>::SetUp() {
  IterationUpdaterAngularFixedUpdaterTest<DimensionWrapper>::SetUp();
  bart::testing::DealiiTestDomain<DimensionWrapper::value>::SetUpDealii();

  auto mock_lhs_ptr = std::make_unique<system::terms::BilinearTermMock>();
  mock_lhs_obs_ptr_ = mock_lhs_ptr.get();

  matrix_ptr_ = std::make_shared<system::MPISparseMatrix>();
  matrix_ptr_->reinit(this->matrix_1);

  test_system_.left_hand_side_ptr_ = std::move(mock_lhs_ptr);
}

TYPED_TEST_SUITE(IterationUpdaterAngularFixedUpdaterDomainTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularFixedUpdaterDomainTest,
    UpdateFixedNullptrMPI) {
  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(this->index_))
      .WillOnce(Return(nullptr));
  EXPECT_ANY_THROW(this->test_updater_ptr_->UpdateFixedTerms(
      this->test_system_,
      this->group_number_,
      this->angle_index_));
}

TYPED_TEST(IterationUpdaterAngularFixedUpdaterDomainTest,
    UpdateFixedMPI) {
  using QuadraturePointType = quadrature::QuadraturePointI<this->dim>;

  // Get three random values to stamp
  auto double_vector = test_helpers::RandomVector(3, 1, 10);
  double sum = std::accumulate(double_vector.begin(), double_vector.end(), 0);
  // Set the value of our expected result
  this->matrix_1 = 0;
  this->StampMatrix(this->matrix_1, sum);

  // Matrix we'll stamp. We'll set it to a value to make sure it gets zero'd
  this->StampMatrix(*this->matrix_ptr_, 100);
  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(this->index_))
      .WillOnce(Return(this->matrix_ptr_));

  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;
  quadrature::QuadraturePointIndex quad_index(this->angle_index_);
  system::EnergyGroup group_number(this->group_number_);

  EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(quad_index))
      .WillOnce(Return(quadrature_point_ptr_));

  auto collision_function = [&](system::MPISparseMatrix& to_stamp) {
    int value = double_vector.at(0);
    this->StampMatrix(to_stamp, value);
  };
  auto streaming_function = [&](system::MPISparseMatrix& to_stamp) {
    int value = double_vector.at(1);
    this->StampMatrix(to_stamp, value);
  };
  auto boundary_function = [&](system::MPISparseMatrix& to_stamp) {
    int value = double_vector.at(2);
    this->StampMatrix(to_stamp, value);
  };
  EXPECT_CALL(*this->stamper_ptr_, StampCollisionTerm(Ref(*this->matrix_ptr_),
                                                      group_number))
      .WillOnce(WithArg<0>(Invoke(collision_function)));
  EXPECT_CALL(*this->stamper_ptr_, StampStreamingTerm(Ref(*this->matrix_ptr_),
                                                      quadrature_point_ptr_,
                                                      group_number))
      .WillOnce(WithArg<0>(Invoke(streaming_function)));
  EXPECT_CALL(*this->stamper_ptr_,
              StampBoundaryBilinearTerm(Ref(*this->matrix_ptr_),
                                        quadrature_point_ptr_,
                                        group_number))
      .WillOnce(WithArg<0>(Invoke(boundary_function)));
  this->test_updater_ptr_->UpdateFixedTerms(this->test_system_,
                                            this->group_number_,
                                            this->angle_index_);
  EXPECT_TRUE(bart::test_helpers::CompareMPIMatrices(*this->matrix_ptr_,
                                                this->matrix_1));
}

} // namespace
