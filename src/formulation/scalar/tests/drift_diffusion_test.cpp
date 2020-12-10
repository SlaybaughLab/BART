#include "formulation/scalar/drift_diffusion.hpp"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "calculator/drift_diffusion/tests/drift_diffusion_vector_calculator_mock.hpp"
#include "data/cross_sections.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;
using test_helpers::AreEqual;

using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::_;


template <typename DimensionWrapper>
class DriftDiffusionFormulationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using CellPtr = typename domain::CellPtr<dim>;
  using CrossSections = data::CrossSections;
  using DriftDiffusionCalculator = NiceMock<typename calculator::drift_diffusion::DriftDiffusionVectorCalculatorMock<dim>>;
  using FiniteElement = NiceMock<typename domain::finite_element::FiniteElementMock<dim>>;
  using Material = NiceMock<btest::MockMaterial>;

  DriftDiffusionFormulationTest() : dof_handler_(triangulation_), finite_element_(1) {};

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
  const int material_id_{ 1 };
  dealii::FullMatrix<double> expected_result_q_0_, expected_result_q_1_;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto DriftDiffusionFormulationTest<DimensionWrapper>::SetUp() -> void {
  drift_diffusion_calculator_mock_ptr_ = std::make_shared<DriftDiffusionCalculator>();
  finite_element_mock_ptr_ = std::make_shared<FiniteElement>();
  Material mock_material;

  ON_CALL(*finite_element_mock_ptr_, dofs_per_cell()).WillByDefault(Return(dofs_per_cell_));
  ON_CALL(*finite_element_mock_ptr_, n_cell_quad_pts()).WillByDefault(Return(cell_quadrature_points_));

  auto ReturnShapeGradient = [&](const int i, const int q) {
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
  };

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    ON_CALL(*finite_element_mock_ptr_, Jacobian(q)).WillByDefault(Return(3 * (q + 1)));
    for (int i = 0; i < dofs_per_cell_; ++i) {
      ON_CALL(*finite_element_mock_ptr_, ShapeValue(i,q)).WillByDefault(Return(1 + i + q));
      ON_CALL(*finite_element_mock_ptr_, ShapeGradient(i,q)).WillByDefault(ReturnShapeGradient);
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

  ON_CALL(*drift_diffusion_calculator_mock_ptr_, DriftDiffusion(_, _, _, _, _)).WillByDefault(Return(drift_diffusion));

  std::unordered_map<int, std::vector<double>> sigma_t{{material_id_, {1.0, 2.0}}};
  std::unordered_map<int, std::vector<double>> diffusion_coef{{material_id_, {0.5, 1.0}}};
  ON_CALL(mock_material, GetSigT()).WillByDefault(Return(sigma_t));
  ON_CALL(mock_material, GetDiffusionCoef()).WillByDefault(Return(sigma_t));
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

  std::array<double, 16> expected_result_q_0_values{0.5, 1.5, 2.5, 3.5,
                                             1.0, 3.0, 5.0, 7.0,
                                             1.5, 4.5, 7.5, 10.5,
                                             2.0, 6.0, 10.0, 14.0};
  std::array<double, 16> expected_result_q_1_values{3.0, 5.0, 7.0, 9.0,
                                             4.5, 7.5, 10.5, 13.5,
                                             6.0, 10.0, 14.0, 18.0,
                                             7.5, 12.5, 17.5, 22.5};
  expected_result_q_0_ = dealii::FullMatrix<double>(2, 2, expected_result_q_0_values.begin());
  expected_result_q_1_ = dealii::FullMatrix<double>(2, 2, expected_result_q_1_values.begin());
}

TYPED_TEST_SUITE(DriftDiffusionFormulationTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionFormulationTest, ConstructorAndDependencyGetters) {
  formulation::scalar::DriftDiffusion<this->dim> drift_diffusion(
      this->finite_element_mock_ptr_,
      this->cross_sections_ptr_,
      this->drift_diffusion_calculator_mock_ptr_);
  EXPECT_NE(drift_diffusion.finite_element_ptr(), nullptr);
  EXPECT_NE(drift_diffusion.cross_sections_ptr(), nullptr);
  EXPECT_NE(drift_diffusion.drift_diffusion_calculator_ptr(), nullptr);
}

TYPED_TEST(DriftDiffusionFormulationTest, ConstructorBadDependencies) {
  const int n_dependencies{ 3 };
  for (int i = 0; i < 3; ++i) {
    EXPECT_ANY_THROW({
      [[maybe_unused]] formulation::scalar::DriftDiffusion<this->dim> drift_diffusion(
          i == 0 ? this->finite_element_mock_ptr_ : nullptr,
          i == 1 ? this->cross_sections_ptr_ : nullptr,
          i == 2 ? this->drift_diffusion_calculator_mock_ptr_ : nullptr);
    });
  }
}



} // namespace
