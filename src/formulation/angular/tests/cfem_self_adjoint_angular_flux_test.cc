#include "formulation/angular/cfem_self_adjoint_angular_flux.h"

#include "data/cross_sections.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Return;

/* Tests for CFEM Self Adjoint Angular Flux formulation class. This class is
 * responsible for filling cell matrices for the SAAF angular formulation.
 *
 * Initial conditions: mock dependencies default call values set; dealii
 * test domain set up; cell pointer initialized to a locally owned cell with
 * material id set to material_id_;
 *
 * Dependencies are all shared pointers so these also function as observation
 * pointers.
 *
 */
template <typename DimensionWrapper>
class FormulationAngularCFEMSelfAdjointAngularFluxTest :
    public ::testing::Test ,
    bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  typename dealii::DoFHandler<dim>::active_cell_iterator cell_ptr_;
  using FiniteElementType = typename domain::finite_element::FiniteElementMock<dim>;
  using QuadratureSetType = typename quadrature::QuadratureSetMock<dim>;
  using MaterialType = btest::MockMaterial;

  // Mock dependencies and supporting objects
  std::shared_ptr<FiniteElementType> mock_finite_element_ptr_;
  std::shared_ptr<QuadratureSetType> mock_quadrature_set_ptr_;
  std::shared_ptr<data::CrossSections> cross_section_ptr_;
  NiceMock<MaterialType> mock_material_;

  // Test parameters
  const int material_id_ = 1;

  void SetUp() override;
};

// SETUP =======================================================================

/* SetUp initial test conditions.
 *
 * 1. Set default values for all mock objects.
 * 2. Set cell pointer to locally owned cell.
 */
template <typename DimensionWrapper>
void FormulationAngularCFEMSelfAdjointAngularFluxTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();

  // Make mocks and set up
  mock_finite_element_ptr_ = std::make_shared<NiceMock<FiniteElementType>>();
  mock_quadrature_set_ptr_ = std::make_shared<NiceMock<QuadratureSetType>>();

  // Set default return values
  ON_CALL(*mock_finite_element_ptr_, dofs_per_cell())
      .WillByDefault(Return(2));
  ON_CALL(*mock_finite_element_ptr_, n_cell_quad_pts())
      .WillByDefault(Return(2));
  ON_CALL(*mock_finite_element_ptr_, n_face_quad_pts())
      .WillByDefault(Return(2));

  // Instantiate cross-section object
  cross_section_ptr_ = std::make_shared<data::CrossSections>(mock_material_);

  // Find an active, locally owned cell and set the material ID
  for (auto cell = this->dof_handler_.begin_active();
       cell != this->dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned()) {
      cell_ptr_ = cell;
      cell_ptr_->set_material_id(material_id_);
      break;
    }
  }
}

TYPED_TEST_CASE(FormulationAngularCFEMSelfAdjointAngularFluxTest,
                bart::testing::AllDimensions);

// TESTS =======================================================================

// CONSTRUCTOR AND CONST ACCESSORS

// Constructor should query appropriate values from finite element object
TYPED_TEST(FormulationAngularCFEMSelfAdjointAngularFluxTest, Constructor) {
  constexpr int dim = this->dim;

  EXPECT_CALL(*this->mock_finite_element_ptr_, dofs_per_cell())
      .WillOnce(::testing::DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, n_cell_quad_pts())
      .WillOnce(::testing::DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, n_face_quad_pts())
      .WillOnce(::testing::DoDefault());

  formulation::angular::CFEMSelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);
}

// Getters should return observation pointers to dependencies
TYPED_TEST(FormulationAngularCFEMSelfAdjointAngularFluxTest,
    DependencyPointers) {
  constexpr int dim = this->dim;

  formulation::angular::CFEMSelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  EXPECT_EQ(test_saaf.finite_element_ptr(), this->mock_finite_element_ptr_.get());
  EXPECT_EQ(test_saaf.cross_sections_ptr(), this->cross_section_ptr_.get());
  EXPECT_EQ(test_saaf.quadrature_set_ptr(), this->mock_quadrature_set_ptr_.get());
}

// FUNCTION TESTS: Initialize
// Initialize should calculate the correct matrices
TYPED_TEST(FormulationAngularCFEMSelfAdjointAngularFluxTest, InitializeValues) {

}

// Initialize should throw an error if cell_ptr is invalid
TYPED_TEST(FormulationAngularCFEMSelfAdjointAngularFluxTest,
    InitializeBadCellPtr) {
  constexpr int dim = this->dim;

  formulation::angular::CFEMSelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::CellPtr<dim> invalid_cell_ptr;
  EXPECT_ANY_THROW(test_saaf.Initialize(invalid_cell_ptr));
}



} // namespace