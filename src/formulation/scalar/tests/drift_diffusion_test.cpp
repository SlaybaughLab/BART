#include "formulation/scalar/drift_diffusion.hpp"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "quadrature/calculators/tests/angular_flux_integrator_mock.hpp"
#include "calculator/drift_diffusion/tests/drift_diffusion_vector_calculator_mock.hpp"
#include "data/cross_sections.h"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;
using test_helpers::AreEqual;

using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::_;
using ::testing::Ref, ::testing::AtLeast;

template <typename DimensionWrapper>
class DriftDiffusionFormulationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using CellPtr = typename domain::CellPtr<dim>;
  using CrossSections = data::CrossSections;
  using AngularFluxIntegrator = NiceMock<quadrature::calculators::AngularFluxIntegratorMock>;
  using DriftDiffusionCalculator = NiceMock<typename calculator::drift_diffusion::DriftDiffusionVectorCalculatorMock<dim>>;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusion<dim>;
  using FiniteElement = NiceMock<typename domain::finite_element::FiniteElementMock<dim>>;
  using Material = NiceMock<btest::MockMaterial>;

  DriftDiffusionFormulationTest() : dof_handler_(triangulation_), finite_element_(1) {};

  // test object
  std::unique_ptr<DriftDiffusionFormulation> test_formulation_{ nullptr };

  std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_mock_ptr_{ nullptr };
  std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_mock_ptr_{ nullptr };
  std::shared_ptr<FiniteElement> finite_element_mock_ptr_{ nullptr };
  std::shared_ptr<CrossSections> cross_sections_ptr_{ nullptr };

  // minimum number of dealii objects to get a valid cell pointer
  CellPtr cell_ptr_;
  dealii::Triangulation<dim> triangulation_;
  dealii::DoFHandler<dim> dof_handler_;
  dealii::FE_Q<dim> finite_element_;

  // Test paraameters
  const int dofs_per_cell_{ 2 };
  const int cell_quadrature_points_{ 2 };
  const int face_quadrature_points_{ 2 };
  const int material_id_{ 1 };
  const int energy_group_{ 0 };
  dealii::FullMatrix<double> expected_result_;
  std::unordered_map<int, std::vector<double>> diffusion_coef_{{material_id_, {0.5, 1.0}}};

  auto SetUp() -> void override;
  [[nodiscard]] auto GetShapeGradient(const int i, const int q) const -> dealii::Tensor<1, dim>;
};


template <typename DimensionWrapper>
auto DriftDiffusionFormulationTest<DimensionWrapper>::GetShapeGradient(const int i, const int q) const
-> dealii::Tensor<1, dim> {
  dealii::Tensor<1, dim> gradient_tensor;
  const double val{ (i + q + 0.5) };
  switch (dim) {
    case 3:
      gradient_tensor[2] = 3 * val; [[fallthrough]];
    case 2:
      gradient_tensor[1] = 2 * val; [[fallthrough]];
    case 1:
      gradient_tensor[0] = val;
      break;
  }
  return gradient_tensor;
}

template <typename DimensionWrapper>
auto DriftDiffusionFormulationTest<DimensionWrapper>::SetUp() -> void {
  angular_flux_integrator_mock_ptr_ = std::make_shared<AngularFluxIntegrator>();
  drift_diffusion_calculator_mock_ptr_ = std::make_shared<DriftDiffusionCalculator>();
  finite_element_mock_ptr_ = std::make_shared<FiniteElement>();
  Material mock_material;

  ON_CALL(*finite_element_mock_ptr_, dofs_per_cell()).WillByDefault(Return(dofs_per_cell_));
  ON_CALL(*finite_element_mock_ptr_, n_cell_quad_pts()).WillByDefault(Return(cell_quadrature_points_));
  ON_CALL(*finite_element_mock_ptr_, n_face_quad_pts()).WillByDefault(Return(face_quadrature_points_));

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    ON_CALL(*finite_element_mock_ptr_, Jacobian(q)).WillByDefault(Return(3 * (q + 1)));
    ON_CALL(*finite_element_mock_ptr_, FaceJacobian(q)).WillByDefault(Return(3 * (q + 1)));
    for (int i = 0; i < dofs_per_cell_; ++i) {
      auto shape_gradient{ GetShapeGradient(i, q) };
      ON_CALL(*finite_element_mock_ptr_, ShapeValue(i,q)).WillByDefault(Return(1 + i + q));
      ON_CALL(*finite_element_mock_ptr_, FaceShapeValue(i,q)).WillByDefault(Return(1 + i + q));
      ON_CALL(*finite_element_mock_ptr_, ShapeGradient(i,q)).WillByDefault(Return(shape_gradient));
    }
  }

  dealii::Tensor<1, dim> drift_diffusion;
  switch (dim) {
    case 3:
      drift_diffusion[2] = 1000; [[fallthrough]];
    case 2:
      drift_diffusion[1] = 100; [[fallthrough]];
    case 1:
      drift_diffusion[0] = 10;
      break;
  }

  ON_CALL(*drift_diffusion_calculator_mock_ptr_, DriftDiffusionVector(_, _, _, _)).WillByDefault(Return(drift_diffusion));

  ON_CALL(mock_material, GetDiffusionCoef()).WillByDefault(Return(diffusion_coef_));
  cross_sections_ptr_ = std::make_shared<CrossSections>(mock_material);

  dealii::GridGenerator::hyper_cube(triangulation_, 0, 1);
  dof_handler_.distribute_dofs(finite_element_);
  for (auto cell = dof_handler_.begin_active(); cell != dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned()) {
      cell_ptr_ = cell;
      cell_ptr_->set_material_id(material_id_);
      break;
    }
  }

  std::array<double, 16> expected_result_q_0_values{0.5, 1.5,
                                                    1.0, 3.0};
  std::array<double, 16> expected_result_q_1_values{3.0, 5.0,
                                                    4.5, 7.5};
  auto expected_result_q_0_ = dealii::FullMatrix<double>(2, 2, expected_result_q_0_values.begin());
  expected_result_q_0_ *= 3;
  auto expected_result_q_1_ = dealii::FullMatrix<double>(2, 2, expected_result_q_1_values.begin());
  expected_result_q_1_ *= 6;
  expected_result_ = expected_result_q_0_;
  expected_result_.add(1, expected_result_q_1_);
  if (dim == 1) {
    expected_result_ *= 10;
  } else if (dim == 2) {
    expected_result_ *= 210;
  } else {
    expected_result_ *= 3210;
  }

  test_formulation_ = std::move(std::make_unique<DriftDiffusionFormulation>(finite_element_mock_ptr_,
                                                                            cross_sections_ptr_,
                                                                            drift_diffusion_calculator_mock_ptr_,
                                                                            angular_flux_integrator_mock_ptr_));
}

TYPED_TEST_SUITE(DriftDiffusionFormulationTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionFormulationTest, ConstructorAndDependencyGetters) {
  EXPECT_CALL(*this->finite_element_mock_ptr_, dofs_per_cell()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_cell_quad_pts()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_face_quad_pts()).WillOnce(DoDefault());
  formulation::scalar::DriftDiffusion<this->dim> drift_diffusion(
      this->finite_element_mock_ptr_,
      this->cross_sections_ptr_,
      this->drift_diffusion_calculator_mock_ptr_,
      this->angular_flux_integrator_mock_ptr_);
  EXPECT_NE(drift_diffusion.finite_element_ptr(), nullptr);
  EXPECT_NE(drift_diffusion.cross_sections_ptr(), nullptr);
  EXPECT_NE(drift_diffusion.drift_diffusion_calculator_ptr(), nullptr);
  EXPECT_NE(drift_diffusion.angular_flux_integrator_ptr(), nullptr);
}

TYPED_TEST(DriftDiffusionFormulationTest, ConstructorBadDependencies) {
  const int n_dependencies{ 4 };
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
      [[maybe_unused]] formulation::scalar::DriftDiffusion<this->dim> drift_diffusion(
          i == 0 ? this->finite_element_mock_ptr_ : nullptr,
          i == 1 ? this->cross_sections_ptr_ : nullptr,
          i == 2 ? this->drift_diffusion_calculator_mock_ptr_ : nullptr,
          i == 3 ? this->angular_flux_integrator_mock_ptr_ : nullptr);
    });
  }
}

TYPED_TEST(DriftDiffusionFormulationTest, FillDriftDiffusion) {
  constexpr int dim = this->dim;
  const int n_dofs(this->dof_handler_.n_dofs());
  auto& finite_element_mock = *this->finite_element_mock_ptr_;
  auto& drift_diffusion_calculator_mock = *this->drift_diffusion_calculator_mock_ptr_;
  std::array<dealii::Vector<double>, dim> current_vectors_at_dofs;
  const dealii::Vector<double> scalar_flux_at_dofs(this->dof_handler_.n_dofs());
  std::vector<double> integrated_angular_flux_at_q{ test_helpers::RandomVector(this->cell_quadrature_points_, 1, 100) };
  std::vector<double> scalar_flux_at_q{ test_helpers::RandomVector(this->cell_quadrature_points_, 1, 100) };
  const double diffusion_coeff{ this->diffusion_coef_.at(this->material_id_).at(this->energy_group_) };

  std::array<std::vector<double>, dim> current_at_q;
  for (auto& current_component : current_at_q)
    current_component = test_helpers::RandomVector(this->cell_quadrature_points_, 1, 100);
  for (int i = 0; i < dim; ++i) {
    current_vectors_at_dofs.at(i) = dealii::Vector<double>(this->dof_handler_.n_dofs());
    auto random_current{ test_helpers::RandomVector(n_dofs, -100, 100) };
    for (int j = 0; j < n_dofs; ++j) {
      current_vectors_at_dofs.at(i)[j] = random_current.at(j);
    }
    EXPECT_CALL(finite_element_mock, ValueAtQuadrature(Ref(current_vectors_at_dofs.at(i))))
        .WillOnce(Return(current_at_q.at(i)));
  }

  EXPECT_CALL(finite_element_mock, ValueAtQuadrature(scalar_flux_at_dofs))
      .WillOnce(Return(scalar_flux_at_q));
  EXPECT_CALL(finite_element_mock, SetCell(this->cell_ptr_));
  for (int q = 0; q < this->cell_quadrature_points_; ++q) {
    EXPECT_CALL(finite_element_mock, Jacobian(q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
    for (int i = 0; i < this->dofs_per_cell_; ++i) {
      EXPECT_CALL(finite_element_mock, ShapeGradient(i, q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
      EXPECT_CALL(finite_element_mock, ShapeValue(i, q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
      dealii::Tensor<1, dim> current;
      for (int j = 0; j < dim; ++j)
        current[j] = current_at_q.at(j).at(q);

      EXPECT_CALL(drift_diffusion_calculator_mock, DriftDiffusionVector(scalar_flux_at_q.at(q),
                                                                        current,
                                                                        this->GetShapeGradient(i, q),
                                                                        diffusion_coeff))
          .Times(AtLeast(1))
          .WillRepeatedly(DoDefault());
    }
  }
  dealii::FullMatrix<double> cell_matrix(this->dofs_per_cell_, this->dofs_per_cell_);
  cell_matrix = 0;
  this->test_formulation_->FillCellDriftDiffusionTerm(cell_matrix,
                                                      this->cell_ptr_,
                                                      system::EnergyGroup(this->energy_group_),
                                                      scalar_flux_at_dofs,
                                                      current_vectors_at_dofs);
  EXPECT_TRUE(AreEqual(this->expected_result_, cell_matrix));
}

TYPED_TEST(DriftDiffusionFormulationTest, FillCellDriftDiffusionBadMatrixSize) {
  std::vector<int> bad_sizes{this->dofs_per_cell_ - 1, this->dofs_per_cell_ + 1};
  const std::array<dealii::Vector<double>, this->dim> current_vectors_at_dofs;
  const dealii::Vector<double> scalar_flux_at_dofs(this->dof_handler_.n_dofs());

  for (const auto bad_size : bad_sizes) {
    dealii::FullMatrix<double> matrix_bad_rows(bad_size, this->dofs_per_cell_);
    dealii::FullMatrix<double> matrix_bad_cols(this->dofs_per_cell_, bad_size);
    for (auto& matrix : std::array<dealii::FullMatrix<double>, 2>{matrix_bad_rows, matrix_bad_cols}) {
      EXPECT_ANY_THROW({
        this->test_formulation_->FillCellDriftDiffusionTerm(matrix,
                                                            this->cell_ptr_,
                                                            system::EnergyGroup(this->energy_group_),
                                                            scalar_flux_at_dofs,
                                                            current_vectors_at_dofs);
      });
    }
  }

}

TYPED_TEST(DriftDiffusionFormulationTest, FillCellBoundaryTerm) {
  constexpr int dim = this->dim;
  using VectorMap = typename formulation::scalar::DriftDiffusion<this->dim>::VectorMap;
  const int n_dofs(this->dof_handler_.n_dofs());
  const domain::FaceIndex face_index{test_helpers::RandomInt(0, 4)};
  auto& finite_element_mock = *this->finite_element_mock_ptr_;

  // Normal values
  dealii::Tensor<1, dim> normal_tensor;
  dealii::Vector<double> normal_vector(dim);

  for (int dir = 0; dir < dim; ++dir) {
    const double random_value{ test_helpers::RandomDouble(-100, 100) };
    normal_tensor[dir] = random_value;
    normal_vector[dir] = random_value;
  }

  // Boundary factor vectors
  std::vector<double> directional_current{ test_helpers::RandomVector(n_dofs, 1, 100) };
  std::vector<double> directional_flux{ test_helpers::RandomVector(n_dofs, 1, 100) };
  dealii::Vector<double> boundary_factor_at_global_dofs(n_dofs);
  const std::vector<double> boundary_factor_at_quadrature{ 10, 20 };

  const int zero_value_position{ test_helpers::RandomInt(0, n_dofs) };
  directional_flux.at(zero_value_position) = 0;

  for (int i = 0; i < n_dofs; ++i)
    boundary_factor_at_global_dofs[i] = directional_current.at(i)/directional_flux.at(i);
  boundary_factor_at_global_dofs[zero_value_position] = 0;

  // Expected results
  const std::array<double, 4> expected_results_values{ 510, 780, 780, 1200 };
  const dealii::FullMatrix<double> expected_results(2, 2, expected_results_values.begin());
  const formulation::BoundaryType reflective_boundary{ formulation::BoundaryType::kVacuum };
  const VectorMap group_angular_flux_map;


  EXPECT_CALL(finite_element_mock, SetFace(this->cell_ptr_, face_index));
  EXPECT_CALL(finite_element_mock, ValueAtFaceQuadrature(::testing::ContainerEq(boundary_factor_at_global_dofs)))
      .WillOnce(Return(boundary_factor_at_quadrature));
  EXPECT_CALL(finite_element_mock, FaceNormal()).WillOnce(Return(normal_tensor));
  EXPECT_CALL(*this->angular_flux_integrator_mock_ptr_, DirectionalCurrent(group_angular_flux_map, normal_vector))
      .WillOnce(Return(directional_current));
  EXPECT_CALL(*this->angular_flux_integrator_mock_ptr_, DirectionalFlux(group_angular_flux_map, normal_vector))
      .WillOnce(Return(directional_flux));

  for (int q = 0; q < this->cell_quadrature_points_; ++q) {
    EXPECT_CALL(finite_element_mock, FaceJacobian(q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
    for (int i = 0; i < this->dofs_per_cell_; ++i) {
      EXPECT_CALL(finite_element_mock, FaceShapeValue(i, q)).Times(AtLeast(2)).WillRepeatedly(DoDefault());
    }
  }
  dealii::FullMatrix<double> cell_matrix(this->dofs_per_cell_, this->dofs_per_cell_);
  cell_matrix = 0;
  this->test_formulation_->FillCellBoundaryTerm(cell_matrix,
                                                this->cell_ptr_,
                                                face_index,
                                                reflective_boundary,
                                                group_angular_flux_map);
  EXPECT_TRUE(AreEqual(expected_results, cell_matrix));
}

TYPED_TEST(DriftDiffusionFormulationTest, FillCellBoundaryTermReflective) {
  using VectorMap = typename formulation::scalar::DriftDiffusion<this->dim>::VectorMap;
  const int n_dofs(this->dof_handler_.n_dofs());
  const domain::FaceIndex face_index{test_helpers::RandomInt(0, 4)};
  dealii::Vector<double> boundary_factor_at_global_dofs(n_dofs);
  const std::array<double, 4> expected_results_values{ 0, 0, 0, 0 };
  const dealii::FullMatrix<double> expected_results(2, 2, expected_results_values.begin());
  const formulation::BoundaryType reflective_boundary{ formulation::BoundaryType::kReflective };

  dealii::FullMatrix<double> cell_matrix(this->dofs_per_cell_, this->dofs_per_cell_);
  cell_matrix = 0;
  this->test_formulation_->FillCellBoundaryTerm(cell_matrix,
                                                this->cell_ptr_,
                                                face_index,
                                                reflective_boundary,
                                                VectorMap());
  EXPECT_TRUE(AreEqual(expected_results, cell_matrix));
}

TYPED_TEST(DriftDiffusionFormulationTest, FillCellBoundaryTermBadMatrixSize) {
  using VectorMap = typename formulation::scalar::DriftDiffusion<this->dim>::VectorMap;
  std::vector<int> bad_sizes{this->dofs_per_cell_ - 1, this->dofs_per_cell_ + 1};
  dealii::Vector<double> boundary_factor_at_global_dofs(this->dof_handler_.n_dofs());
  const domain::FaceIndex face_index{test_helpers::RandomInt(0, 4)};
  const formulation::BoundaryType reflective_boundary{ formulation::BoundaryType::kVacuum };
  for (const auto bad_size : bad_sizes) {
    dealii::FullMatrix<double> matrix_bad_rows(bad_size, this->dofs_per_cell_);
    dealii::FullMatrix<double> matrix_bad_cols(this->dofs_per_cell_, bad_size);
    for (auto& matrix : std::array<dealii::FullMatrix<double>, 2>{matrix_bad_rows, matrix_bad_cols}) {
      EXPECT_ANY_THROW({
                         this->test_formulation_->FillCellBoundaryTerm(matrix,
                                                                       this->cell_ptr_,
                                                                       face_index,
                                                                       reflective_boundary,
                                                                       VectorMap());
                       });
    }
  }
}

} // namespace
