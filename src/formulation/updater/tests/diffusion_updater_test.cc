#include "formulation/updater/diffusion_updater.h"

#include "formulation/scalar/tests/diffusion_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace {

using namespace bart;

using ::testing::DoDefault, ::testing::_, ::testing::Ref, ::testing::Invoke,
::testing::WithArg;

template <typename DimensionWrapper>
class FormulationUpdaterDiffusionTest :
    public bart::formulation::updater::test_helpers::UpdaterTests<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;
  using Boundary = problem::Boundary;

  std::unique_ptr<UpdaterType> test_updater_ptr_;

  FormulationType* formulation_obs_ptr_;
  StamperType* stamper_obs_ptr_;

  std::unordered_set<Boundary> reflective_boundaries{Boundary::kXMin, Boundary::kYMax};

  void SetUp() override;
};

TYPED_TEST_SUITE(FormulationUpdaterDiffusionTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void FormulationUpdaterDiffusionTest<DimensionWrapper>::SetUp() {
  bart::formulation::updater::test_helpers::UpdaterTests<dim>::SetUp();
  auto formulation_ptr = std::make_unique<FormulationType>();
  formulation_obs_ptr_ = formulation_ptr.get();
  auto stamper_ptr = this->MakeStamper();
  stamper_obs_ptr_ = stamper_ptr.get();

  test_updater_ptr_ = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                    std::move(stamper_ptr),
                                                    reflective_boundaries);
}

// ===== CONSTRUCTOR TESTS =====================================================

TYPED_TEST(FormulationUpdaterDiffusionTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_NO_THROW(
      {
        test_updater_ptr = std::make_unique<UpdaterType>(
            std::move(formulation_ptr), std::move(stamper_ptr));
      }
  );
  ASSERT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  ASSERT_NE(test_updater_ptr->stamper_ptr(), nullptr);
}

TYPED_TEST(FormulationUpdaterDiffusionTest, ConstructorReflective) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unordered_set<problem::Boundary> reflective_boundaries{problem::Boundary::kXMin};
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_NO_THROW(
      {
        test_updater_ptr = std::make_unique<UpdaterType>(
            std::move(formulation_ptr), std::move(stamper_ptr),
            reflective_boundaries);
      }
  );
  ASSERT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  ASSERT_NE(test_updater_ptr->stamper_ptr(), nullptr);
  EXPECT_EQ(reflective_boundaries, test_updater_ptr->reflective_boundaries());
}

TYPED_TEST(FormulationUpdaterDiffusionTest, ConstructorBadDependencies) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
      nullptr, std::move(stamper_ptr)); });
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
        std::move(formulation_ptr), nullptr); });
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
      nullptr, nullptr); });
}

// ===== UpdateFixedTerms() TESTS ==============================================

TYPED_TEST(FormulationUpdaterDiffusionTest, UpdateFixedTermTest) {
  system::EnergyGroup group_number(this->group_number);
  quadrature::QuadraturePointIndex angle_index(this->angle_index);
  bart::system::Index scalar_index{this->group_number, 0};

  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(scalar_index))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetFixedTermPtr(scalar_index))
      .WillOnce(DoDefault());

  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellStreamingTerm(_, cell, this->group_number));
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellCollisionTerm(_, cell, this->group_number));
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellFixedSource(_, cell, this->group_number));
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          problem::Boundary boundary_id = static_cast<problem::Boundary>(
              cell->face(face)->boundary_id());

          using BoundaryType = typename formulation::scalar::DiffusionI<this->dim>::BoundaryType;
          BoundaryType boundary_type = BoundaryType::kVacuum;
          if (this->reflective_boundaries.count(boundary_id) == 1)
            boundary_type = BoundaryType::kReflective;

          EXPECT_CALL(*this->formulation_obs_ptr_,
                      FillBoundaryTerm(_, cell, face, boundary_type));
        }
      }
    }
  }

  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampMatrix(Ref(*this->matrix_to_stamp), _))
      .Times(2)
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampBoundaryMatrix(Ref(*this->matrix_to_stamp), _))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampVector(Ref(*this->vector_to_stamp), _))
      .WillOnce(DoDefault());

  this->test_updater_ptr_->UpdateFixedTerms(this->test_system_, group_number, angle_index);
  EXPECT_TRUE(test_helpers::CompareMPIMatrices(this->expected_result,
                                               *this->matrix_to_stamp));
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

// ====== UpdateScatteringSource TESTS =========================================

TYPED_TEST(FormulationUpdaterDiffusionTest, UpdateScatteringSourceTest) {
  system::EnergyGroup group_number(this->group_number);
  quadrature::QuadraturePointIndex angle_index(this->angle_index);
  bart::system::Index scalar_index{this->group_number, 0};

  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetVariableTermPtr(
      scalar_index, system::terms::VariableLinearTerms::kScatteringSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(_,_))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(DoDefault());
  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_, FillCellScatteringSource(
        _, cell, group_number.get(), Ref(this->current_iteration_moments_)))
        .WillOnce(DoDefault());
  }

  this->test_updater_ptr_->UpdateScatteringSource(this->test_system_, group_number, angle_index);
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

// ===== UpdateFissionSource TEST ==============================================
TYPED_TEST(FormulationUpdaterDiffusionTest, UpdateFissionSourceTest) {
  system::EnergyGroup group_number(this->group_number);
  quadrature::QuadraturePointIndex angle_index(this->angle_index);
  bart::system::Index scalar_index{this->group_number, 0};

  const double k_effective = bart::test_helpers::RandomDouble(0, 1.5);
  this->test_system_.k_effective = k_effective;

  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetVariableTermPtr(
      scalar_index, system::terms::VariableLinearTerms::kFissionSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(_,_))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(DoDefault());
  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_, FillCellFissionSource(
        _, cell, group_number.get(), k_effective,
        Ref(this->current_iteration_moments_.at({group_number.get(), 0, 0})),
        Ref(this->current_iteration_moments_)))
        .WillOnce(DoDefault());
  }

  this->test_updater_ptr_->UpdateFissionSource(this->test_system_, group_number, angle_index);
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

// ===== UpdateFixedSource TEST ================================================
TYPED_TEST(FormulationUpdaterDiffusionTest, UpdateFixedSourceTest) {
  system::EnergyGroup group_number(this->group_number);
  quadrature::QuadraturePointIndex angle_index(this->angle_index);
  bart::system::Index scalar_index{this->group_number, 0};

  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetFixedTermPtr(scalar_index))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(_,_))
      .WillOnce(DoDefault());

  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
        FillCellFixedSource(_, cell, group_number.get()))
        .WillOnce(DoDefault());
  }

  this->test_updater_ptr_->UpdateFixedSource(this->test_system_, group_number, angle_index);
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

} // namespace