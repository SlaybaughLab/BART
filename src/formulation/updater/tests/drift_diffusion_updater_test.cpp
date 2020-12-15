#include "formulation/updater/drift_diffusion_updater.hpp"

#include <deal.II/lac/vector.h>

#include "quadrature/calculators/tests/drift_diffusion_integrated_flux_mock.hpp"
#include "formulation/scalar/tests/diffusion_mock.h"
#include "formulation/scalar/tests/drift_diffusion_mock.hpp"
#include "formulation/tests/stamper_mock.h"
#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"
#include "system/solution/solution_types.h"
#include "system/moments/tests/spherical_harmonic_mock.h"

namespace  {

using namespace bart;
template <int dim>
using UpdaterTest = bart::formulation::updater::test_helpers::UpdaterTests<dim>;


using ::testing::ContainerEq, ::testing::DoDefault, ::testing::_, ::testing::Ref, ::testing::Return, ::testing::ReturnRef;

template <typename DimensionWrapper>
class FormulationUpdaterDriftDiffusionTest : public UpdaterTest<DimensionWrapper::value> {
 public:
  static constexpr int dim{ DimensionWrapper::value };
  using AngularFluxStorageMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using DiffusionFormulation = formulation::scalar::DiffusionMock<dim>;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusionMock<dim>;
  using IntegratedFluxCalculator = quadrature::calculators::DriftDiffusionIntegratedFluxMock;
  using Stamper = formulation::StamperMock<dim>;
  using Boundary = problem::Boundary;
  using Vector = dealii::Vector<double>;
  using GroupAngularFluxMap = std::map<quadrature::QuadraturePointIndex, std::shared_ptr<Vector>>;
  using HighOrderMoments = system::moments::SphericalHarmonicMock;

  // Test object
  using TestUpdater = formulation::updater::DriftDiffusionUpdater<dim>;
  std::unique_ptr<TestUpdater> test_updater_ptr_{ nullptr };

  // Supporting objects
  AngularFluxStorageMap angular_flux_storage_map_{}; // All angular fluxes
  GroupAngularFluxMap group_angular_flux_{}; // Angular fluxes for a single group (the test group)
  dealii::Vector<double> integrated_angular_flux_;
  dealii::Vector<double> group_scalar_flux_;

  // Mock observation pointers
  DiffusionFormulation* diffusion_formulation_obs_ptr_{ nullptr };
  DriftDiffusionFormulation* drift_diffusion_formulation_obs_ptr_{ nullptr };
  IntegratedFluxCalculator* integrated_flux_calculator_obs_ptr_{ nullptr };
  Stamper* stamper_obs_ptr_{ nullptr };
  std::shared_ptr<HighOrderMoments> high_order_moments_ptr_{ nullptr };

  // Test parameters (each should be a different value to ensure the correct values are being passed
  const Vector::size_type angular_flux_size{ static_cast<unsigned int>(test_helpers::RandomInt(5, 10)) };

  std::unordered_set<Boundary> reflective_boundaries{Boundary::kXMin, Boundary::kYMax};

  void SetUp() override;
};

template<typename DimensionWrapper>
void FormulationUpdaterDriftDiffusionTest<DimensionWrapper>::SetUp() {
  UpdaterTest<dim>::SetUp();
  auto diffusion_formulation_ptr = std::make_unique<DiffusionFormulation>();
  diffusion_formulation_obs_ptr_ = diffusion_formulation_ptr.get();
  auto drift_diffusion_formulation_ptr = std::make_unique<DriftDiffusionFormulation>();
  drift_diffusion_formulation_obs_ptr_ = drift_diffusion_formulation_ptr.get();
  auto integrated_flux_calculator_ptr = std::make_unique<IntegratedFluxCalculator>();
  integrated_flux_calculator_obs_ptr_ = integrated_flux_calculator_ptr.get();
  auto stamper_ptr = std::shared_ptr<Stamper>(this->MakeStamper());
  stamper_obs_ptr_ = stamper_ptr.get();
  high_order_moments_ptr_ = std::make_shared<HighOrderMoments>();

  integrated_angular_flux_ = dealii::Vector<double>(angular_flux_size);
  group_scalar_flux_ = dealii::Vector<double>(angular_flux_size);


  using Index = system::SolutionIndex;
  for (int group = 0; group < this->total_groups; ++group) {
    for (int angle = 0; angle < this->total_angles; ++angle) {
      auto angular_flux = std::make_shared<Vector>(angular_flux_size);
      Index solution_index{ system::EnergyGroup(group), system::AngleIdx(angle) };
      angular_flux_storage_map_.insert({solution_index, angular_flux});
      if (group == this->group_number)
        group_angular_flux_.insert({quadrature::QuadraturePointIndex(angle), angular_flux});
    }
  }
  test_updater_ptr_ = std::make_unique<TestUpdater>(std::move(diffusion_formulation_ptr),
                                                    std::move(drift_diffusion_formulation_ptr),
                                                    stamper_ptr,
                                                    std::move(integrated_flux_calculator_ptr),
                                                    high_order_moments_ptr_,
                                                    angular_flux_storage_map_,
                                                    reflective_boundaries);
}

TYPED_TEST_SUITE(FormulationUpdaterDriftDiffusionTest, bart::testing::AllDimensions);

// Constructor should not throw, dependencies should have been saved and getters return non-null ptrs
TYPED_TEST(FormulationUpdaterDriftDiffusionTest, ConstructorDependencyGetters) {
  constexpr int dim{ this->dim };
  using DiffusionFormulation = formulation::scalar::DiffusionMock<dim>;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusionMock<dim>;
  using IntegratedFluxCalculator = quadrature::calculators::DriftDiffusionIntegratedFluxMock;
  using Updater = formulation::updater::DriftDiffusionUpdater<dim>;
  using Stamper = formulation::StamperMock<dim>;

  std::unique_ptr<Updater> test_updater{ nullptr };
  EXPECT_NO_THROW({
    test_updater = std::make_unique<Updater>(std::make_unique<DiffusionFormulation>(),
                                             std::make_unique<DriftDiffusionFormulation>(),
                                             std::make_shared<Stamper>(),
                                             std::make_unique<IntegratedFluxCalculator>(),
                                             this->high_order_moments_ptr_,
                                             this->angular_flux_storage_map_);
  });
  ASSERT_NE(test_updater->formulation_ptr(), nullptr);
  ASSERT_NE(test_updater->drift_diffusion_formulation_ptr(), nullptr);
  ASSERT_NE(test_updater->stamper_ptr(), nullptr);
  EXPECT_EQ(test_updater->high_order_moments(), this->high_order_moments_ptr_.get());
  ASSERT_THAT(test_updater->angular_flux_storage_map(), ContainerEq(this->angular_flux_storage_map_));
}

TYPED_TEST(FormulationUpdaterDriftDiffusionTest, ConstructorBadDependencies) {
  constexpr int dim{ this->dim };
  using DiffusionFormulation = formulation::scalar::DiffusionMock<dim>;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusionMock<dim>;
  using IntegratedFluxCalculator = quadrature::calculators::DriftDiffusionIntegratedFluxMock;
  using Updater = formulation::updater::DriftDiffusionUpdater<dim>;
  using Stamper = formulation::StamperMock<dim>;
  using HighOrderMoments = system::moments::SphericalHarmonicMock;


  std::unique_ptr<Updater> test_updater{ nullptr };
  const int n_dependencies_to_check{ 5 };
  for (int i = 0; i < n_dependencies_to_check; ++i) {
    EXPECT_ANY_THROW({
      test_updater = std::make_unique<Updater>(i == 0 ? nullptr : std::make_unique<DiffusionFormulation>(),
                                               i == 1 ? nullptr : std::make_unique<DriftDiffusionFormulation>(),
                                               i == 2 ? nullptr : std::make_unique<Stamper>(),
                                               i == 3 ? nullptr : std::make_unique<IntegratedFluxCalculator>(),
                                               i == 4 ? nullptr : std::make_shared<HighOrderMoments>(),
                                               this->angular_flux_storage_map_);
                     });
  }
}

TYPED_TEST(FormulationUpdaterDriftDiffusionTest, UpdateFixedTermTest) {
  system::EnergyGroup group_number(this->group_number);
  quadrature::QuadraturePointIndex angle_index(this->angle_index);
  bart::system::Index scalar_index{this->group_number, 0};

  EXPECT_CALL(*this->mock_lhs_obs_ptr_, GetFixedTermPtr(scalar_index)).WillOnce(DoDefault());
  EXPECT_CALL(*this->mock_rhs_obs_ptr_, GetFixedTermPtr(scalar_index)).WillOnce(DoDefault());
  EXPECT_CALL(*this->integrated_flux_calculator_obs_ptr_, Integrate(ContainerEq(this->group_angular_flux_)))
      .WillOnce(Return(this->integrated_angular_flux_));
  std::array<int, 3> moment_index{this->group_number, 0, 0};
  EXPECT_CALL(*this->high_order_moments_ptr_, GetMoment(moment_index)).WillOnce(ReturnRef(this->group_scalar_flux_));

  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->diffusion_formulation_obs_ptr_, FillCellStreamingTerm(_, cell, this->group_number));
    EXPECT_CALL(*this->diffusion_formulation_obs_ptr_, FillCellCollisionTerm(_, cell, this->group_number));
    EXPECT_CALL(*this->diffusion_formulation_obs_ptr_, FillCellFixedSource(_, cell, this->group_number));
    EXPECT_CALL(*this->drift_diffusion_formulation_obs_ptr_,
        FillCellDriftDiffusionTerm(_,
                                   cell,
                                   system::EnergyGroup(this->group_number),
                                   this->group_scalar_flux_,
                                   this->integrated_angular_flux_));
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          problem::Boundary boundary_id = static_cast<problem::Boundary>(cell->face(face)->boundary_id());

          using BoundaryType = typename formulation::scalar::DiffusionI<this->dim>::BoundaryType;
          BoundaryType boundary_type = BoundaryType::kVacuum;
          if (this->reflective_boundaries.count(boundary_id) == 1)
            boundary_type = BoundaryType::kReflective;

          EXPECT_CALL(*this->diffusion_formulation_obs_ptr_, FillBoundaryTerm(_, cell, face, boundary_type));
        }
      }
    }
  }

  EXPECT_CALL(*this->stamper_obs_ptr_, StampMatrix(Ref(*this->matrix_to_stamp), _))
      .Times(3)
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_, StampBoundaryMatrix(Ref(*this->matrix_to_stamp), _)).WillOnce(DoDefault());
  EXPECT_CALL(*this->stamper_obs_ptr_, StampVector(Ref(*this->vector_to_stamp), _)).WillOnce(DoDefault());

  this->test_updater_ptr_->UpdateFixedTerms(this->test_system_, group_number, angle_index);
  EXPECT_TRUE(test_helpers::AreEqual(this->expected_result, *this->matrix_to_stamp));
  EXPECT_TRUE(test_helpers::AreEqual(this->expected_vector_result, *this->vector_to_stamp));
}

} // namespace
