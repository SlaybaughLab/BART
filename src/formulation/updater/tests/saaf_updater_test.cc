#include "formulation/updater/saaf_updater.h"

#include "quadrature/tests/quadrature_set_mock.h"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "system/terms/tests/bilinear_term_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.h"

namespace {

using namespace bart;

using ::testing::Return, ::testing::Ref, ::testing::Invoke, ::testing::_,
::testing::A, ::testing::WithArg;

template <typename DimensionWrapper>
class FormulationUpdaterSAAFTest :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  // Test object
  std::unique_ptr<UpdaterType> test_updater_ptr;

  // Pointers to mocks
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  FormulationType* formulation_obs_ptr_;
  StamperType* stamper_obs_ptr_;

  // Other test objects and parameters
  system::System test_system_;
  std::shared_ptr<bart::system::MPISparseMatrix> matrix_to_stamp;
  bart::system::terms::BilinearTermMock* mock_lhs_obs_ptr_;
  system::MPISparseMatrix& expected_result = this->matrix_2;
  const int group_number = test_helpers::RandomDouble(0, 10);
  const int angle_index = test_helpers::RandomDouble(0, 10);
  const system::Index index{group_number, angle_index};
  static domain::CellPtr<dim> static_cell_ptr;

  static void EvaluateFunction(std::function<void(formulation::FullMatrix&,
                                                  const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateBoundaryFunction(std::function<void(formulation::FullMatrix&,
                                                          const domain::FaceIndex,
                                                          const domain::CellPtr<dim>&)> stamp_function);
  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationUpdaterSAAFTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  auto formulation_ptr = std::make_unique<FormulationType>();
  formulation_obs_ptr_ = formulation_ptr.get();
  auto stamper_ptr = std::make_unique<StamperType>();
  stamper_obs_ptr_ = stamper_ptr.get();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();
  test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                   std::move(stamper_ptr),
                                                   quadrature_set_ptr_);

  auto mock_lhs_ptr = std::make_unique<system::terms::BilinearTermMock>();
  mock_lhs_obs_ptr_ = mock_lhs_ptr.get();
  matrix_to_stamp = std::make_shared<system::MPISparseMatrix>();
  matrix_to_stamp->reinit(this->matrix_1);

  test_system_.left_hand_side_ptr_ = std::move(mock_lhs_ptr);
}

template <typename DimensionWrapper>
void FormulationUpdaterSAAFTest<DimensionWrapper>::EvaluateFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  stamp_function(to_stamp, cell_ptr);
}

template <typename DimensionWrapper>
void FormulationUpdaterSAAFTest<DimensionWrapper>::EvaluateBoundaryFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  domain::FaceIndex index(0);
  stamp_function(to_stamp, index, cell_ptr);
}

TYPED_TEST_SUITE(FormulationUpdaterSAAFTest, bart::testing::AllDimensions);

// ====== Constructor tests ====================================================

TYPED_TEST(FormulationUpdaterSAAFTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;

  EXPECT_NO_THROW({
    test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                     std::move(stamper_ptr),
                                                     quadrature_set_ptr);
  });
  EXPECT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->stamper_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->quadrature_set_ptr(), nullptr);
}

TYPED_TEST(FormulationUpdaterSAAFTest, ConstructorBadDepdendencies) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  for (bool formulation_good : {true, false}) {
    for (bool stamper_good : {true, false}) {
      for (bool quadrature_set_good : {true, false}) {
        auto formulation_ptr = formulation_good ?
                               std::make_unique<FormulationType>() : nullptr;
        auto stamper_ptr = stamper_good ?
                           std::make_unique<StamperType>() : nullptr;
        auto quadrature_set_ptr = quadrature_set_good ?
            std::make_shared<QuadratureSetType>() : nullptr;
        if (!formulation_good || !stamper_good || !quadrature_set_good) {
          std::unique_ptr<UpdaterType> test_updater_ptr;
          EXPECT_ANY_THROW({
            test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                std::move(stamper_ptr), quadrature_set_ptr);
          });
        }
      }
    }
  }
}



// ===== Update Fixed Terms Tests ==============================================

TYPED_TEST(FormulationUpdaterSAAFTest, UpdateFixedTermsTest) {
  constexpr int dim = this->dim;
  using QuadraturePointType = quadrature::QuadraturePointI<dim>;
  // Stamp the matrix we are going to stamp with a random value, it should be
  // reset to zero before being passed to the stamper.
  this->StampMatrix(*this->matrix_to_stamp, test_helpers::RandomDouble(-10, 10));
  this->expected_result = 0;
  quadrature::QuadraturePointIndex quad_index(this->angle_index);
  system::EnergyGroup group_number(this->group_number);
  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;

  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(this->index))
      .WillOnce(Return(this->matrix_to_stamp));
  EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(quad_index))
      .WillOnce(Return(quadrature_point_ptr_));

  EXPECT_CALL(*this->formulation_obs_ptr_,
              FillCellStreamingTerm(_, _, quadrature_point_ptr_, group_number));
  EXPECT_CALL(*this->formulation_obs_ptr_,
              FillCellCollisionTerm(_, _, group_number));
  EXPECT_CALL(*this->formulation_obs_ptr_,
              FillBoundaryBilinearTerm(_, _, _, quadrature_point_ptr_, group_number));

  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampMatrix(Ref(*this->matrix_to_stamp),_))
      .Times(2)
      .WillRepeatedly(WithArg<1>(Invoke(this->EvaluateFunction)));

  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampBoundaryMatrix(Ref(*this->matrix_to_stamp),_))
      .WillOnce(WithArg<1>(Invoke(this->EvaluateBoundaryFunction)));

  this->test_updater_ptr->UpdateFixedTerms(this->test_system_, group_number,
                                           quad_index);
  EXPECT_TRUE(test_helpers::CompareMPIMatrices(this->expected_result,
                                               *this->matrix_to_stamp));
}

} // namespace
