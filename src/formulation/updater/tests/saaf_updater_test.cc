#include "formulation/updater/saaf_updater.h"

#include "quadrature/tests/quadrature_set_mock.h"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"
#include "system/solution/solution_types.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"

namespace {

using namespace bart;

using ::testing::Return, ::testing::Ref, ::testing::Invoke, ::testing::_,
::testing::A, ::testing::WithArg, ::testing::DoDefault, ::testing::ReturnRef,
::testing::NiceMock;

template <typename DimensionWrapper>
class FormulationUpdaterSAAFTest :
    public bart::formulation::updater::test_helpers::UpdaterTests<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using AngularSolutionType = NiceMock<system::solution::MPIGroupAngularSolutionMock>;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;
  using Boundary = problem::Boundary;

  // Test object
  std::unique_ptr<UpdaterType> test_updater_ptr;

  // Pointers to mocks
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  FormulationType* formulation_obs_ptr_;
  StamperType* stamper_obs_ptr_;

  std::unordered_set<Boundary> reflective_boundaries{Boundary::kXMin,
                                                     Boundary::kYMax,
                                                     Boundary::kZMin};
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map_;
  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationUpdaterSAAFTest<DimensionWrapper>::SetUp() {
  bart::formulation::updater::test_helpers::UpdaterTests<dim>::SetUp();
  auto formulation_ptr = std::make_unique<FormulationType>();
  formulation_obs_ptr_ = formulation_ptr.get();
  auto stamper_ptr = this->MakeStamper();
  stamper_obs_ptr_ = stamper_ptr.get();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  auto angular_solution_ptr = std::make_shared<AngularSolutionType>();
  ON_CALL(*angular_solution_ptr, total_angles())
      .WillByDefault(Return(this->total_angles));

  angular_solution_ptr_map_.insert({system::EnergyGroup(this->group_number),
                                    angular_solution_ptr});

  ON_CALL(*quadrature_set_ptr_, size())
      .WillByDefault(Return(this->total_angles));

  test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                   std::move(stamper_ptr),
                                                   quadrature_set_ptr_,
                                                   angular_solution_ptr_map_,
                                                   reflective_boundaries);
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
  EXPECT_EQ(test_updater_ptr->reflective_boundaries().size(), 0);
}

TYPED_TEST(FormulationUpdaterSAAFTest, ConstructorReflective) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using AngularSolutionType = bart::system::solution::MPIGroupAngularSolutionMock;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map;

  auto angular_solution_ptr = std::make_shared<AngularSolutionType>();

  EXPECT_CALL(*quadrature_set_ptr, size()).WillOnce(Return(this->total_angles));

  for (int group = 0; group < this->total_groups; ++group) {
    auto angular_solution_ptr = std::make_shared<AngularSolutionType>();
    EXPECT_CALL(*angular_solution_ptr, total_angles())
        .WillOnce(Return(this->total_angles));
    angular_solution_ptr_map.insert(
        {bart::system::EnergyGroup(group), angular_solution_ptr});
  }

  EXPECT_NO_THROW({
    test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                     std::move(stamper_ptr),
                                                     quadrature_set_ptr,
                                                     angular_solution_ptr_map,
                                                     this->reflective_boundaries);
                  });
  EXPECT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->stamper_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->quadrature_set_ptr(), nullptr);
  EXPECT_EQ(test_updater_ptr->angular_solution_ptr_map(),
            angular_solution_ptr_map);
  EXPECT_EQ(test_updater_ptr->reflective_boundaries(), this->reflective_boundaries);
}

TYPED_TEST(FormulationUpdaterSAAFTest, ConstructorReflectiveBadAngleMatch) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using AngularSolutionType = bart::system::solution::MPIGroupAngularSolutionMock;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map;

  EXPECT_CALL(*quadrature_set_ptr, size()).WillOnce(Return(this->total_angles));

  for (int group = 0; group < this->total_groups; ++group) {
    auto angular_solution_ptr = std::make_shared<AngularSolutionType>();
    if (group != this->group_number) {
      EXPECT_CALL(*angular_solution_ptr, total_angles())
          .Times(::testing::AtMost(1))
          .WillRepeatedly(Return(this->total_angles));
    } else {
      EXPECT_CALL(*angular_solution_ptr, total_angles())
          .WillOnce(Return(this->total_angles + 1));
    }
    angular_solution_ptr_map.insert(
        {bart::system::EnergyGroup(group), angular_solution_ptr});
  }

  EXPECT_ANY_THROW({
    test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                     std::move(stamper_ptr),
                                                     quadrature_set_ptr,
                                                     angular_solution_ptr_map,
                                                     this->reflective_boundaries);
                  });
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

  quadrature::QuadraturePointIndex quad_index(this->angle_index);
  system::EnergyGroup group_number(this->group_number);
  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;

  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(this->index))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(quad_index))
      .WillOnce(Return(quadrature_point_ptr_));

  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellStreamingTerm(_,
                                      cell,
                                      quadrature_point_ptr_,
                                      group_number));
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellCollisionTerm(_, cell, group_number));
    int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
    if (cell->at_boundary()) {
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          EXPECT_CALL(*this->formulation_obs_ptr_,
                      FillBoundaryBilinearTerm(_, cell, domain::FaceIndex(face),
                                               quadrature_point_ptr_,
                                               group_number));
        }
      }
    }
  }

  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampMatrix(Ref(*this->matrix_to_stamp),_))
      .Times(2)
      .WillRepeatedly(DoDefault());

  EXPECT_CALL(*this->stamper_obs_ptr_,
      StampBoundaryMatrix(Ref(*this->matrix_to_stamp),_))
      .WillOnce(DoDefault());

  this->test_updater_ptr->UpdateFixedTerms(this->test_system_, group_number,
                                           quad_index);
  EXPECT_TRUE(test_helpers::CompareMPIMatrices(this->expected_result,
                                               *this->matrix_to_stamp));
}

TYPED_TEST(FormulationUpdaterSAAFTest, UpdateScatteringSourceTest) {
  constexpr int dim = this->dim;
  using QuadraturePointType = quadrature::QuadraturePointI<dim>;

  quadrature::QuadraturePointIndex quad_index(this->angle_index);
  system::EnergyGroup group_number(this->group_number);
  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;

  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetVariableTermPtr(
      this->index,
      system::terms::VariableLinearTerms::kScatteringSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(quad_index))
      .WillOnce(Return(quadrature_point_ptr_));
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(_,_))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(DoDefault());
  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_, FillCellScatteringSourceTerm(
        _, cell, quadrature_point_ptr_, group_number,
        Ref(this->current_iteration_moments_.at({group_number.get(), 0, 0})),
        Ref(this->current_iteration_moments_)));
  }

  this->test_updater_ptr->UpdateScatteringSource(this->test_system_,
                                                 group_number, quad_index);
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

TYPED_TEST(FormulationUpdaterSAAFTest, UpdateFissionSourceTest) {
  constexpr int dim = this->dim;
  using QuadraturePointType = quadrature::QuadraturePointI<dim>;

  quadrature::QuadraturePointIndex quad_index(this->angle_index);
  system::EnergyGroup group_number(this->group_number);
  std::shared_ptr<QuadraturePointType> quadrature_point_ptr_;

  const double k_effective = 1.045;
  this->test_system_.k_effective = k_effective;

  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetVariableTermPtr(
      this->index,
      system::terms::VariableLinearTerms::kFissionSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(quad_index))
      .WillOnce(Return(quadrature_point_ptr_));
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(_,_))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->current_moments_obs_ptr_, moments())
      .WillOnce(DoDefault());
  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_, FillCellFissionSourceTerm(
        _, cell, quadrature_point_ptr_, group_number, k_effective,
        Ref(this->current_iteration_moments_.at({group_number.get(), 0, 0})),
        Ref(this->current_iteration_moments_)));
  }

  this->test_updater_ptr->UpdateFissionSource(this->test_system_,
                                              group_number, quad_index);
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->expected_vector_result,
                                              *this->vector_to_stamp));
}

} // namespace
