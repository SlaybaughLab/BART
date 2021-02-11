#include "formulation/angular/self_adjoint_angular_flux.h"

#include <deal.II/base/tensor.h>

#include "data/cross_sections.h"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "material/tests/material_mock.hpp"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "quadrature/utility/quadrature_utilities.h"
#include "system/system_types.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_assertions.hpp"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;
using ::testing::ByRef;
using ::testing::_;
using ::testing::A;

namespace test_helpers = bart::test_helpers;
using test_helpers::AreEqual;

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
class FormulationAngularSelfAdjointAngularFluxTest : public ::testing::Test ,
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

  // Other test objects
  std::set<std::shared_ptr<quadrature::QuadraturePointI<dim>>,
           quadrature::utility::quadrature_point_compare<dim>> quadrature_set_;
  std::set<int> quadrature_point_indices_ = {};

  // Test parameters
  const int material_id_ = 1, non_fissile_material_id_ = 0;
  const double k_effective_ = 1.107;
  // This factor is used to get fission transfer matrices from sigma_s_per_ster_
  const double fission_test_factor_ = 2.3465;
  // Cross-sections
  const std::unordered_map<int, std::vector<double>> sigma_t_{{material_id_, {1.0, 2.0}}};
  const std::unordered_map<int, std::vector<double>> inv_sigma_t_{{material_id_, {1.0, 0.5}},
                                                                  {non_fissile_material_id_, {1.0, 0.5}}};
  const std::unordered_map<int, std::vector<double>> q_per_ster_{{non_fissile_material_id_, {1.0, 2.0}}};
  const std::unordered_map<int, formulation::FullMatrix> sigma_s_per_ster_{
    {material_id_,
     {2, 2, std::array<double, 4>{0.25, 0.5, 0.75, 1.0}.begin()}}};
  const std::unordered_map<int, formulation::FullMatrix> fission_xfer_per_ster_{
      {material_id_,
       {2,2, std::array<double, 4>{0.25*fission_test_factor_,
                                   0.75*fission_test_factor_,
                                   0.5*fission_test_factor_,
                                   fission_test_factor_}.begin()}}};
  system::moments::MomentVector group_0_moment_, group_1_moment_;
  system::moments::MomentsMap out_group_moments_;
  const std::vector<double> group_0_moment_values_{0.75, 0.75};
  const std::vector<double> group_1_moment_values_{1.0, 1.0};

  void SetUp() override;
};

// SETUP =======================================================================

/* SetUp initial test conditions.
 *
 * 1. Set default values for all mock objects.
 * 2. Set cell pointer to locally owned cell.
 */
template <typename DimensionWrapper>
void FormulationAngularSelfAdjointAngularFluxTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();

  group_0_moment_.reinit(2);
  group_1_moment_.reinit(2);

  // Make mocks and set up
  mock_finite_element_ptr_ = std::make_shared<NiceMock<FiniteElementType>>();
  mock_quadrature_set_ptr_ = std::make_shared<NiceMock<QuadratureSetType>>();

  // Set default return values
  ON_CALL(*mock_finite_element_ptr_, dofs_per_cell()).WillByDefault(Return(2));
  ON_CALL(*mock_finite_element_ptr_, n_cell_quad_pts()).WillByDefault(Return(2));
  ON_CALL(*mock_finite_element_ptr_, n_face_quad_pts()).WillByDefault(Return(2));
  // Set default return values for shape functions and Jacobian
  for (int quad_pt_idx = 0; quad_pt_idx < 2; ++quad_pt_idx) {
    int jacobian = 3*(quad_pt_idx + 1);
    ON_CALL(*mock_finite_element_ptr_, Jacobian(quad_pt_idx)).WillByDefault(Return(jacobian));
    ON_CALL(*mock_finite_element_ptr_, FaceJacobian(quad_pt_idx)).WillByDefault(Return(jacobian));
    for (int dof_idx = 0; dof_idx < 2; ++dof_idx) {
      int entry = quad_pt_idx + 1 + 10*(dof_idx + 1);
      ON_CALL(*mock_finite_element_ptr_, ShapeValue(dof_idx, quad_pt_idx)).WillByDefault(Return(entry));
      ON_CALL(*mock_finite_element_ptr_, FaceShapeValue(dof_idx, quad_pt_idx)).WillByDefault(Return(entry));
      dealii::Tensor<1, dim> gradient_entry;
      for (int i = 0; i < dim; ++i)
        gradient_entry[i] = (entry);
      ON_CALL(*mock_finite_element_ptr_, ShapeGradient(dof_idx, quad_pt_idx)).WillByDefault(Return(gradient_entry));
    }
  }

  // Set up mock quadrature points for quadrature set
  for (int n_angle = 0; n_angle < 2; ++n_angle) {
    auto new_quadrature_point =
        std::make_shared<NiceMock<quadrature::QuadraturePointMock<dim>>>();
    std::array<double, dim> position;
    position.fill(n_angle + 1);
    dealii::Tensor<1, dim> tensor_position;
    for (int i = 0; i < dim; ++i)
      tensor_position[i] = (n_angle + 1);
    ON_CALL(*new_quadrature_point, cartesian_position()).WillByDefault(Return(position));
    ON_CALL(*new_quadrature_point, cartesian_position_tensor()).WillByDefault(Return(tensor_position));
    ON_CALL(*mock_quadrature_set_ptr_,GetQuadraturePoint(quadrature::QuadraturePointIndex(n_angle)))
        .WillByDefault(Return(new_quadrature_point));

    auto return_pair = quadrature_set_.insert(new_quadrature_point);
    quadrature_point_indices_.insert(n_angle);

    ON_CALL(*mock_quadrature_set_ptr_, GetQuadraturePointIndex(*return_pair.first)).WillByDefault(Return(n_angle));
  }

  ON_CALL(*mock_quadrature_set_ptr_, quadrature_point_indices()).WillByDefault(Return(quadrature_point_indices_));
  ON_CALL(*mock_quadrature_set_ptr_, begin()).WillByDefault(Return(quadrature_set_.begin()));
  ON_CALL(*mock_quadrature_set_ptr_, end()).WillByDefault(Return(quadrature_set_.end()));

  // Set up cross-sections
  ON_CALL(mock_material_, GetSigT()).WillByDefault(Return(sigma_t_));
  ON_CALL(mock_material_, GetInvSigT()).WillByDefault(Return(inv_sigma_t_));
  ON_CALL(mock_material_, GetQPerSter()).WillByDefault(Return(q_per_ster_));
  ON_CALL(mock_material_, GetSigSPerSter()).WillByDefault(Return(sigma_s_per_ster_));
  ON_CALL(mock_material_, GetChiNuSigFPerSter()).WillByDefault(Return(fission_xfer_per_ster_));
  ON_CALL(mock_material_, GetFissileIDMap()).WillByDefault(Return(std::unordered_map<int, bool>{
        {material_id_, true},
        {non_fissile_material_id_, false}}));

  // Set up moment values
  for (int i = 0; i < 2; ++i) {
    group_0_moment_(i) = group_0_moment_values_.at(i);
    group_1_moment_(i) = group_1_moment_values_.at(i);
  }

  out_group_moments_[{0,0,0}] = group_0_moment_;
  out_group_moments_[{1,0,0}] = group_1_moment_;

  ON_CALL(*mock_finite_element_ptr_, ValueAtQuadrature(group_0_moment_)).WillByDefault(Return(group_0_moment_values_));
  ON_CALL(*mock_finite_element_ptr_, ValueAtQuadrature(group_1_moment_)).WillByDefault(Return(group_1_moment_values_));
  ON_CALL(*mock_finite_element_ptr_, ValueAtFaceQuadrature(_)).WillByDefault(Return(group_0_moment_values_));


  // Instantiate cross-section object
  cross_section_ptr_ = std::make_shared<data::CrossSections>(mock_material_);

  // Find an active, locally owned cell and set the material ID
  for (auto cell = this->dof_handler_.begin_active(); cell != this->dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned()) {
      cell_ptr_ = cell;
      cell_ptr_->set_material_id(material_id_);
      break;
    }
  }
}

TYPED_TEST_CASE(FormulationAngularSelfAdjointAngularFluxTest, bart::testing::AllDimensions);

// TESTS =======================================================================

// CONSTRUCTOR AND CONST ACCESSORS

// Constructor should query appropriate values from finite element object
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, Constructor) {
  constexpr int dim = this->dim;

  EXPECT_CALL(*this->mock_finite_element_ptr_, dofs_per_cell()).WillOnce(::testing::DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, n_cell_quad_pts()).WillOnce(::testing::DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, n_face_quad_pts()).WillOnce(::testing::DoDefault());

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);
  EXPECT_FALSE(test_saaf.is_initialized());
}

// Getters should return observation pointers to dependencies
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, DependencyPointers) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  EXPECT_EQ(test_saaf.finite_element_ptr(), this->mock_finite_element_ptr_.get());
  EXPECT_EQ(test_saaf.cross_sections_ptr(), this->cross_section_ptr_.get());
  EXPECT_EQ(test_saaf.quadrature_set_ptr(), this->mock_quadrature_set_ptr_.get());
  EXPECT_FALSE(test_saaf.is_initialized());
}

// FUNCTION TESTS: Initialize

// Initialize should calculate the correct matrices for shape squared
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, InitializeValuesShapeSquared) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  // Quadrature point shape-squared matrix expected values based on
  // hand-calculating the values from shape(i,q) = q + 1 + 10*(i + 1)
  formulation::FullMatrix expected_shape_squared_q_0{2,2, std::array<double, 4>{121, 231, 231, 441}.begin()};
  formulation::FullMatrix expected_shape_squared_q_1{2,2, std::array<double, 4>{144, 264, 264, 484}.begin()};

  EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_)).Times(1);
  EXPECT_CALL(*this->mock_finite_element_ptr_, ShapeValue(_,_)).Times(16).WillRepeatedly(DoDefault());

  EXPECT_NO_THROW(test_saaf.Initialize(this->cell_ptr_));
  auto shape_squared = test_saaf.shape_squared();

  ASSERT_EQ(shape_squared.size(), 2);
  EXPECT_TRUE(AreEqual(expected_shape_squared_q_0, shape_squared.at(0)));
  EXPECT_TRUE(AreEqual(expected_shape_squared_q_1, shape_squared.at(1)));
  EXPECT_TRUE(test_saaf.is_initialized());
}

// Initialize should calculate the correct matrices for Omega dot gradient vectors
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, InitializeOmegaDotGradient) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  /* Procedure should be: get all the quadrature point indices, retrieve each
   * quadrature point using indices, get the position tensor and multiply by
   * the gradient shape. This should be repeated for each degree of freedom
   * (i.e. twice). */
  EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_)).Times(1);
  EXPECT_CALL(*this->mock_quadrature_set_ptr_, quadrature_point_indices()).Times(::testing::AtLeast(1))
      .WillRepeatedly(DoDefault());

  for (int index : this->quadrature_point_indices_) {
    EXPECT_CALL(*this->mock_quadrature_set_ptr_,GetQuadraturePoint(quadrature::QuadraturePointIndex(index)))
        .Times(2)
        .WillRepeatedly(DoDefault());
  }

  for (auto quadrature_point_ptr : this->quadrature_set_) {
    auto mock_quadrature_point_ptr =dynamic_cast<quadrature::QuadraturePointMock<dim>*>(quadrature_point_ptr.get());
    EXPECT_CALL(*mock_quadrature_point_ptr, cartesian_position_tensor())
        .Times(4)
        .WillRepeatedly(DoDefault());
  }

  for (int quad_pt_idx = 0; quad_pt_idx < 2; ++quad_pt_idx) {
    for (int dof_idx = 0; dof_idx < 2; ++dof_idx) {
      EXPECT_CALL(*this->mock_finite_element_ptr_,
          ShapeGradient(quad_pt_idx, dof_idx))
          .Times(2)
          .WillRepeatedly(DoDefault());
    }
  }

  /* Expected values: For each quadrature point (with a given omega) there should
   * be two entries, a value for each degree of freedom. */

  std::map<int, std::map<int, std::vector<double>>> omega_dot_gradient{
      {0, {{0, {11*dim, 21*dim}}, {1, {22*dim, 42*dim}}}},
      {1, {{0, {12*dim, 22*dim}}, {1, {24*dim, 44*dim}}}}
  };

  EXPECT_NO_THROW(test_saaf.Initialize(this->cell_ptr_));

  for (int cell_quad_point = 0; cell_quad_point < 2; ++cell_quad_point) {
    for (int angle_index : this->quadrature_point_indices_) {
      std::vector<double> result;
      EXPECT_NO_THROW(result = test_saaf.OmegaDotGradient(cell_quad_point,
                                                           quadrature::QuadraturePointIndex(angle_index)));
      EXPECT_THAT(result, ::testing::ContainerEq(omega_dot_gradient.at(cell_quad_point).at(angle_index)));
    }
  }
  EXPECT_TRUE(test_saaf.is_initialized());
}

/* Initialize should generate the correct squares of the omega dot gradient
 * vectors. Returned as a matrix with entries corresponding to cell quadrature
 * points (i,j). No extra expectations are required, these are covered by
 * the above test. */
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,InitializeOmegaDotGradientSquared) {
  constexpr int dim = this->dim;

  double dim_sq = dim*dim;
  formulation::FullMatrix odgs_00{2,2, std::array<double, 4>{121, 231, 231, 441}.begin()};
  odgs_00 *= dim_sq;
  formulation::FullMatrix odgs_01{2,2, std::array<double, 4>{484, 924, 924, 1764}.begin()};
  odgs_01 *= dim_sq;
  formulation::FullMatrix odgs_10{2,2, std::array<double, 4>{144, 264, 264, 484}.begin()};
  odgs_10 *= dim_sq;
  formulation::FullMatrix odgs_11{2,2, std::array<double, 4>{576, 1056, 1056, 1936}.begin()};
  odgs_11 *= dim_sq;
  std::map<int, std::map<int, formulation::FullMatrix>> omega_dot_gradient_squared = {
      {0, {{0, odgs_00}, {1, odgs_01}}},
      {1, {{0, odgs_10}, {1, odgs_11}}}
      };

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);
  test_saaf.Initialize(this->cell_ptr_);

  for (int cell_quad_point = 0; cell_quad_point < 2; ++cell_quad_point) {
    for (int angle_index : this->quadrature_point_indices_) {
      formulation::FullMatrix result;
      ASSERT_NO_THROW(result = test_saaf.OmegaDotGradientSquared(cell_quad_point,
                                                                 quadrature::QuadraturePointIndex(angle_index)));
      EXPECT_TRUE(AreEqual(omega_dot_gradient_squared.at(cell_quad_point).at(angle_index), result));
    }
  }
  EXPECT_TRUE(test_saaf.is_initialized());
}

// Initialize should throw an error if cell_ptr is invalid
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, InitializeBadCellPtr) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  domain::CellPtr<dim> invalid_cell_ptr;
  EXPECT_ANY_THROW(test_saaf.Initialize(invalid_cell_ptr));
  EXPECT_FALSE(test_saaf.is_initialized());
}
// =============================================================================
// FUNCTION TESTS:
// =============================================================================

// FillBoundaryBilinearTerm
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillBoundaryBilinearTermTestBadCellPtr) {
  const int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);
  auto angle_ptr = *this->quadrature_set_.begin();
  domain::CellPtr<dim> invalid_cell_ptr;

  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({
    test_saaf.FillBoundaryBilinearTerm(cell_matrix, invalid_cell_ptr,
                                       domain::FaceIndex(0), angle_ptr,
                                       system::EnergyGroup(0));
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellBoundaryBilinearTermBadMatrixSize) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix bad_cell_matrix(3,2);
  auto angle_ptr = *this->quadrature_set_.begin();


  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({
    test_saaf.FillBoundaryBilinearTerm(bad_cell_matrix, this->cell_ptr_,
                                       domain::FaceIndex(0), angle_ptr,
                                       system::EnergyGroup(0));
                   });

  formulation::FullMatrix second_bad_cell_matrix(2,3);
  EXPECT_ANY_THROW({
    test_saaf.FillBoundaryBilinearTerm(second_bad_cell_matrix, 
                                       this->cell_ptr_, domain::FaceIndex(0),
                                       angle_ptr, system::EnergyGroup(0));
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellBoundaryTermLessThanZeroTest) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix expected_results(2,2);
  formulation::FullMatrix cell_matrix(2, 2, std::array<double, 4>{2454, 4554, 4554, 8454}.begin());
  expected_results = cell_matrix;

  test_saaf.Initialize(this->cell_ptr_);
  auto angle_ptr = *this->quadrature_set_.begin();
  int face_index = 0;

  dealii::Tensor<1, dim> normal;
  for (int i = 0; i < dim; ++i)
    normal[i] = -1;
  EXPECT_CALL(*this->mock_finite_element_ptr_, SetFace(this->cell_ptr_, domain::FaceIndex(face_index)));
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceNormal()).WillOnce(Return(normal));

  EXPECT_NO_THROW({
  test_saaf.FillBoundaryBilinearTerm(cell_matrix, this->cell_ptr_,
                                     domain::FaceIndex(0), angle_ptr,
                                     system::EnergyGroup(0));
                  });
  EXPECT_TRUE(AreEqual(expected_results, cell_matrix));
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellBoundaryTermTest) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix expected_results(2, 2, std::array<double, 4>{1227, 2277, 2277, 4227}.begin());
  expected_results *= 3*dim;
  formulation::FullMatrix cell_matrix(2,2);
  cell_matrix = 0;

  test_saaf.Initialize(this->cell_ptr_);
  auto angle_ptr = *this->quadrature_set_.begin();
  int face_index = 0;

  dealii::Tensor<1, dim> normal;
  for (int i = 0; i < dim; ++i)
    normal[i] = 3;

  EXPECT_CALL(*this->mock_finite_element_ptr_, SetFace(this->cell_ptr_, domain::FaceIndex(face_index)));
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceNormal()).WillOnce(Return(normal));
  auto mock_angle_ptr = dynamic_cast<quadrature::QuadraturePointMock<dim>*>(angle_ptr.get());
  EXPECT_CALL(*mock_angle_ptr, cartesian_position_tensor()).WillOnce(DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceJacobian(_)).Times(2).WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceShapeValue(_, _)).Times(16).WillRepeatedly(DoDefault());

  test_saaf.FillBoundaryBilinearTerm(cell_matrix, this->cell_ptr_, domain::FaceIndex(0), angle_ptr,
                                     system::EnergyGroup(0));

  EXPECT_TRUE(AreEqual(expected_results, cell_matrix));
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellBoundaryTermTestNotInitialized) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);

  auto angle_ptr = *this->quadrature_set_.begin();
  EXPECT_ANY_THROW({
  test_saaf.FillBoundaryBilinearTerm(cell_matrix, this->cell_ptr_,
                                     domain::FaceIndex(0), angle_ptr,
                                     system::EnergyGroup(0));
                   });
}

// FillReflectiveBoundaryLinearTerm ======================================================

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillReflectiveBoundaryLinearTermTestBadCellPtr) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  domain::CellPtr<dim> invalid_cell_ptr;
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillReflectiveBoundaryLinearTerm(cell_vector, invalid_cell_ptr, domain::FaceIndex(0), angle_ptr,
                                               dealii::Vector<double>{});
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillReflectiveBoundaryLinearTermTestBadVectorSize) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(3);
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
                     test_saaf.FillReflectiveBoundaryLinearTerm(cell_vector,
                                                      this->cell_ptr_,
                                                      domain::FaceIndex(0),
                                                      angle_ptr,
                                                      dealii::Vector<double>{});
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillReflectiveBoundaryLinearTermTestGreaterThanZero) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector expected_results(2);
  formulation::Vector cell_vector(2);
  cell_vector[0] = test_helpers::RandomDouble(1, 1000);
  cell_vector[1] = test_helpers::RandomDouble(1, 1000);
  expected_results = cell_vector;

  test_saaf.Initialize(this->cell_ptr_);
  auto angle_ptr = *this->quadrature_set_.begin();
  int face_index = 0;

  dealii::Tensor<1, dim> normal;

  EXPECT_CALL(*this->mock_finite_element_ptr_, SetFace(this->cell_ptr_, domain::FaceIndex(face_index)));
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceNormal())
      .WillOnce(Return(normal));

  EXPECT_NO_THROW({
                    test_saaf.FillReflectiveBoundaryLinearTerm(cell_vector,
                                                     this->cell_ptr_,
                                                     domain::FaceIndex(0),
                                                     angle_ptr,
                                                     dealii::Vector<double>{});
                  });

  EXPECT_EQ(expected_results, cell_vector);
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillReflectiveBoundaryLinearTermTest) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector expected_results(2);
  expected_results[0] = 78.75 * dim;
  expected_results[1] = 146.25 * dim;
  formulation::Vector cell_vector(2);

  test_saaf.Initialize(this->cell_ptr_);
  auto angle_ptr = *this->quadrature_set_.begin();
  int face_index = 0;

  dealii::Tensor<1, dim> normal;
  for (int i = 0; i < dim; ++i)
    normal[i] = -1;
  EXPECT_CALL(*this->mock_finite_element_ptr_, SetFace(this->cell_ptr_, domain::FaceIndex(face_index)));
  EXPECT_CALL(*this->mock_finite_element_ptr_, FaceNormal()).WillOnce(Return(normal));

  for (const auto f_q : {0, 1}) {
    EXPECT_CALL(*this->mock_finite_element_ptr_, FaceJacobian(f_q)).WillOnce(DoDefault());
    for (const auto dof : {0, 1}) {
      EXPECT_CALL(*this->mock_finite_element_ptr_, FaceShapeValue(dof, f_q)).WillOnce(DoDefault());
    }
  }


  EXPECT_CALL(*this->mock_finite_element_ptr_, ValueAtFaceQuadrature(_)).WillOnce(DoDefault());

  EXPECT_NO_THROW({
                    test_saaf.FillReflectiveBoundaryLinearTerm(cell_vector, this->cell_ptr_, domain::FaceIndex(0),
                                                               angle_ptr, dealii::Vector<double>{});
                  });

  EXPECT_EQ(expected_results, cell_vector);
}

// FillStreamingTerm ===========================================================
TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellStreamingTermTestBadCellPtr) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);
  auto angle_ptr = *this->quadrature_set_.begin();
  domain::CellPtr<dim> invalid_cell_ptr;


  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({
    test_saaf.FillCellStreamingTerm(cell_matrix, invalid_cell_ptr, angle_ptr, system::EnergyGroup(0));
  });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellStreamingTermTestBadMatrixSize) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix bad_cell_matrix(3,2);
  auto angle_ptr = *this->quadrature_set_.begin();


  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({
    test_saaf.FillCellStreamingTerm(bad_cell_matrix, this->cell_ptr_,
                                    angle_ptr, system::EnergyGroup(0));
                   });

  formulation::FullMatrix second_bad_cell_matrix(2,3);
  EXPECT_ANY_THROW({
    test_saaf.FillCellStreamingTerm(second_bad_cell_matrix, 
                                    this->cell_ptr_, angle_ptr,
                                    system::EnergyGroup(0));
                   });
}


TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellStreamingTermTest) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);
  test_saaf.Initialize(this->cell_ptr_);

  std::map<std::pair<int, int>, formulation::FullMatrix> expected_results;
  formulation::FullMatrix expected_result_g0_a0(2, 2, std::array<double, 4>{1227, 2277, 2277, 4227}.begin());
  formulation::FullMatrix expected_result_g1_a0(2, 2, std::array<double, 4>{613.5, 1138.5, 1138.5, 2113.5}.begin());
  formulation::FullMatrix expected_result_g0_a1(2, 2, std::array<double, 4>{4908, 9108, 9108, 16908}.begin());
  formulation::FullMatrix expected_result_g1_a1(2, 2, std::array<double, 4>{2454, 4554, 4554, 8454}.begin());
  expected_results.insert_or_assign({0,0}, expected_result_g0_a0);
  expected_results.insert_or_assign({1,0}, expected_result_g1_a0);
  expected_results.insert_or_assign({0,1}, expected_result_g0_a1);
  expected_results.insert_or_assign({1,1}, expected_result_g1_a1);

  for (auto& result : expected_results)
    result.second *= dim*dim;

  for (int group = 0; group < 2; ++group) {
    for (int angle = 0; angle < 2; ++angle) {
      auto angle_it = this->quadrature_set_.begin();
      if (angle == 1)
        ++angle_it;
      auto angle_ptr = *angle_it;
      EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_));
      EXPECT_CALL(*this->mock_finite_element_ptr_, Jacobian(_)).Times(2).WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_quadrature_set_ptr_, GetQuadraturePointIndex(_)).WillOnce(DoDefault());

      cell_matrix = 0;

      EXPECT_NO_THROW({
        test_saaf.FillCellStreamingTerm(cell_matrix, this->cell_ptr_,
                                        angle_ptr, system::EnergyGroup(group));
                      });
      std::pair<int, int> result_index{group, angle};
      EXPECT_TRUE(AreEqual(expected_results.at(result_index), cell_matrix))
                      << "Failed: group: " << group << " angle: " << angle;
    }
  }
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellStreamingTermTestNotInitialized) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);

  auto angle_ptr = *this->quadrature_set_.begin();
  EXPECT_ANY_THROW({
    test_saaf.FillCellStreamingTerm(cell_matrix, this->cell_ptr_,
                                    angle_ptr, system::EnergyGroup(0));
                   });
}

// FillCollisionTerm ===========================================================

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellCollisionTermTestBadCellPtr) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);
  domain::CellPtr<dim> invalid_cell_ptr;
  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({test_saaf.FillCellCollisionTerm(cell_matrix, 
                                                    invalid_cell_ptr,
                                                    system::EnergyGroup(0));});
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellCollisionTermTestNotInitialized) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix cell_matrix(2,2);
  EXPECT_ANY_THROW({test_saaf.FillCellCollisionTerm(cell_matrix, this->cell_ptr_, system::EnergyGroup(0));});
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellCollisionTermTestBadMatrixSize) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::FullMatrix bad_cell_matrix(3,2), second_bad_cell_matrix(2,3);;
  test_saaf.Initialize(this->cell_ptr_);
  EXPECT_ANY_THROW({test_saaf.FillCellCollisionTerm(bad_cell_matrix, this->cell_ptr_, system::EnergyGroup(0));});
  EXPECT_ANY_THROW({test_saaf.FillCellCollisionTerm(second_bad_cell_matrix, this->cell_ptr_, system::EnergyGroup(0));});
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellCollisionTermTest) {
  formulation::angular::SelfAdjointAngularFlux<this->dim> test_saaf(this->mock_finite_element_ptr_,
                                                                    this->cross_section_ptr_,
                                                                    this->mock_quadrature_set_ptr_);

  test_saaf.Initialize(this->cell_ptr_);

  for (int group = 0; group < 2; ++group) {
    EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_));
    EXPECT_CALL(*this->mock_finite_element_ptr_, Jacobian(_)).Times(2).WillRepeatedly(DoDefault());
    EXPECT_CALL(*this->mock_finite_element_ptr_, ShapeValue(_,_)).Times(16).WillRepeatedly(DoDefault());

    formulation::FullMatrix cell_matrix(2,2), expected_result(2,2,
                                                              std::array<double, 4>{1227, 2277, 2277, 4227}.begin());
    expected_result *= this->sigma_t_.at(this->material_id_).at(group);

    EXPECT_NO_THROW({test_saaf.FillCellCollisionTerm(cell_matrix, this->cell_ptr_,system::EnergyGroup(group));});
    EXPECT_TRUE(AreEqual(expected_result, cell_matrix));
  }
}

// FillCellFixedSourceTerm =====================================================

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellFixedSourceTermBadCell) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  domain::CellPtr<dim> invalid_cell_ptr;
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellFixedSourceTerm(cell_vector, invalid_cell_ptr, angle_ptr, system::EnergyGroup(0));
  });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellFixedSourceTermTestNotInitialized) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  auto angle_ptr = *this->quadrature_set_.begin();

  EXPECT_ANY_THROW({
    test_saaf.FillCellFixedSourceTerm(cell_vector, this->cell_ptr_, angle_ptr, system::EnergyGroup(0));
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellFixedSourceTermBadVectorLength) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector bad_cell_vector(3);
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellFixedSourceTerm(bad_cell_vector, this->cell_ptr_, angle_ptr, system::EnergyGroup(0));
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellFixedSourceTermFissileMaterial) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2), expected_results(2);
  expected_results = 0;
  test_saaf.Initialize(this->cell_ptr_);

  for (int group = 0; group < 2; ++group) {
    for (int angle = 0; angle < 2; ++angle) {
      auto angle_it = this->quadrature_set_.begin();
      if (angle == 1)
        ++angle_it;
      auto angle_ptr = *angle_it;
      cell_vector = 0;

      EXPECT_NO_THROW({
        test_saaf.FillCellFixedSourceTerm(cell_vector, this->cell_ptr_, angle_ptr, system::EnergyGroup(group));
                      });
      std::pair<int, int> result_index{group, angle};
      EXPECT_EQ(expected_results, cell_vector) << "Failed: group: " << group << " angle: " << angle;
    }
  }
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellFixedSourceTerm) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(this->mock_finite_element_ptr_,
                                                              this->cross_section_ptr_,
                                                              this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  test_saaf.Initialize(this->cell_ptr_);
  this->cell_ptr_->set_material_id(this->non_fissile_material_id_);

  double d = this->dim;
  std::map<std::pair<int, int>, formulation::Vector> expected_results;

  formulation::Vector expected_result_g0_a0(2), expected_result_g0_a1(2),
      expected_result_g1_a0(2), expected_result_g1_a1(2);
  expected_result_g0_a0[0] = 105.0 * (d + 1.0);
  expected_result_g0_a0[1] = 195.0*(d + 1.0);
  expected_result_g0_a1[0] = 105.0 + 210.0*d;
  expected_result_g0_a1[1] = 195.0 + 390.0*d;
  expected_result_g1_a0[0] = 210.0 + 105.0*d;
  expected_result_g1_a0[1] = 390.0 + 195.0*d;
  expected_result_g1_a1[0] = 210.0*(d + 1.0);
  expected_result_g1_a1[1] = 390.0*(d + 1.0);

  expected_results.insert_or_assign({0,0}, expected_result_g0_a0);
  expected_results.insert_or_assign({1,0}, expected_result_g1_a0);
  expected_results.insert_or_assign({0,1}, expected_result_g0_a1);
  expected_results.insert_or_assign({1,1}, expected_result_g1_a1);

  for (int group = 0; group < 2; ++group) {
    for (int angle = 0; angle < 2; ++angle) {
      auto angle_it = this->quadrature_set_.begin();
      if (angle == 1)
        ++angle_it;
      auto angle_ptr = *angle_it;
      EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_));
      EXPECT_CALL(*this->mock_finite_element_ptr_, Jacobian(_))
          .Times(2)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_finite_element_ptr_, ShapeValue(_,_))
          .Times(4)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_quadrature_set_ptr_, GetQuadraturePointIndex(_))
          .WillOnce(DoDefault());

      cell_vector = 0;

      //EXPECT_NO_THROW({
        test_saaf.FillCellFixedSourceTerm(cell_vector, this->cell_ptr_, angle_ptr, system::EnergyGroup(group));
//                      });
      std::pair<int, int> result_index{group, angle};
      EXPECT_EQ(expected_results.at(result_index), cell_vector)
                << "Failed: group: " << group << " angle: " << angle;
    }
  }
}

// FillCellScatteringSourceTerm =====================================================

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellScatteringSourceTermBadCell) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  domain::CellPtr<dim> invalid_cell_ptr;
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellScatteringSourceTerm(cell_vector, invalid_cell_ptr,
                                           angle_ptr, system::EnergyGroup(0),
                                           this->group_0_moment_,
                                           this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest, FillCellScatteringSourceTermTestNotInitialized) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  auto angle_ptr = *this->quadrature_set_.begin();

  EXPECT_ANY_THROW({
    test_saaf.FillCellScatteringSourceTerm(cell_vector, this->cell_ptr_,
                                           angle_ptr, system::EnergyGroup(0),
                                           this->group_0_moment_,
                                           this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellScatteringSourceTermBadVectorLength) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector bad_cell_vector(3);
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellScatteringSourceTerm(bad_cell_vector, this->cell_ptr_,
                                           angle_ptr, system::EnergyGroup(0),
                                           this->group_0_moment_,
                                           this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellScatteringSourceTerm) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  test_saaf.Initialize(this->cell_ptr_);

  double d = this->dim;

  std::map<std::pair<int, int>, formulation::Vector> expected_results;

  formulation::Vector expected_result_g0_a0(2), expected_result_g0_a1(2),
      expected_result_g1_a0(2), expected_result_g1_a1(2);
  expected_result_g0_a0[0] = 72.1875 + 72.1875*d;
  expected_result_g0_a0[1] = 134.0625 * (1 + d);
  expected_result_g0_a1[0] = 72.1875 + 144.375*d;
  expected_result_g0_a1[1] = 134.0625 + 268.125*d;
  expected_result_g1_a0[0] = 164.0625 + 82.03125*d;
  expected_result_g1_a0[1] = 304.6875 + 152.34375*d;
  expected_result_g1_a1[0] = 164.0625 * (d + 1.0);
  expected_result_g1_a1[1] = 304.6875 * (d + 1.0);
  expected_results.insert_or_assign({0,0}, expected_result_g0_a0);
  expected_results.insert_or_assign({1,0}, expected_result_g1_a0);
  expected_results.insert_or_assign({0,1}, expected_result_g0_a1);
  expected_results.insert_or_assign({1,1}, expected_result_g1_a1);

  for (int group = 0; group < 2; ++group) {
    for (int angle = 0; angle < 2; ++angle) {
      auto angle_it = this->quadrature_set_.begin();
      if (angle == 1)
        ++angle_it;
      auto angle_ptr = *angle_it;
      EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_));
      EXPECT_CALL(*this->mock_finite_element_ptr_, Jacobian(_))
          .Times(2)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_finite_element_ptr_, ShapeValue(_,_))
          .Times(4)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_quadrature_set_ptr_, GetQuadraturePointIndex(_))
          .WillOnce(DoDefault());

      int out_group = !group;
      std::array<int, 3> out_index{out_group, 0, 0}, in_index{group, 0, 0};
      auto out_group_moments_map = this->out_group_moments_;
      out_group_moments_map.at(in_index) = 0;

      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(this->out_group_moments_.at(out_index)))
          .Times(1)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(this->out_group_moments_.at(in_index)))
          .Times(0);


      system::moments::MomentVector& in_group_moment = this->group_0_moment_;
      if (group == 1)
        in_group_moment = this->group_1_moment_;

      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(in_group_moment))
          .Times(1)
          .WillRepeatedly(DoDefault());

      cell_vector = 0;

      EXPECT_NO_THROW({
        test_saaf.FillCellScatteringSourceTerm(cell_vector, this->cell_ptr_,
                                               angle_ptr,
                                               system::EnergyGroup(group),
                                               in_group_moment,
                                               out_group_moments_map);
                                        });

      std::pair<int, int> result_index{group, angle};
      EXPECT_EQ(expected_results.at(result_index), cell_vector)
                << "Failed: group: " << group << " angle: " << angle;
    }
  }
}

// FillCellFissionSourceTerm =====================================================

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellFissionSourceTermBadCell) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  domain::CellPtr<dim> invalid_cell_ptr;
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellFissionSourceTerm(cell_vector, invalid_cell_ptr,
                                        angle_ptr, system::EnergyGroup(0),
                                        this->k_effective_,
                                        this->group_0_moment_,
                                        this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellFissionSourceTermTestNotInitialized) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  auto angle_ptr = *this->quadrature_set_.begin();

  EXPECT_ANY_THROW({
    test_saaf.FillCellFissionSourceTerm(cell_vector, this->cell_ptr_,
                                        angle_ptr, system::EnergyGroup(0),
                                        this->k_effective_,
                                        this->group_0_moment_,
                                        this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellFissionSourceTermBadVectorLength) {
  constexpr int dim = this->dim;
  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector bad_cell_vector(3);
  auto angle_ptr = *this->quadrature_set_.begin();
  test_saaf.Initialize(this->cell_ptr_);

  EXPECT_ANY_THROW({
    test_saaf.FillCellFissionSourceTerm(bad_cell_vector, this->cell_ptr_,
                                        angle_ptr, system::EnergyGroup(0),
                                        this->k_effective_,
                                        this->group_0_moment_,
                                        this->out_group_moments_);
                   });
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellFissionSourceTermNonFissileMaterial) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  test_saaf.Initialize(this->cell_ptr_);
  this->cell_ptr_->set_material_id(this->non_fissile_material_id_);
  auto angle_ptr = *this->quadrature_set_.begin();
  const int group = 0;
  system::moments::MomentVector& in_group_moment = this->group_0_moment_;
  auto out_group_moments_map = this->out_group_moments_;

  EXPECT_NO_THROW({
        test_saaf.FillCellFissionSourceTerm(cell_vector, this->cell_ptr_,
                                            angle_ptr,
                                            system::EnergyGroup(group),
                                            this->k_effective_,
                                            in_group_moment,
                                            out_group_moments_map);
                      });
  EXPECT_TRUE(cell_vector.l2_norm() < 1e-10);
}

TYPED_TEST(FormulationAngularSelfAdjointAngularFluxTest,
           FillCellFissionSourceTerm) {
  constexpr int dim = this->dim;

  formulation::angular::SelfAdjointAngularFlux<dim> test_saaf(
      this->mock_finite_element_ptr_,
      this->cross_section_ptr_,
      this->mock_quadrature_set_ptr_);

  formulation::Vector cell_vector(2);
  test_saaf.Initialize(this->cell_ptr_);

  double d = this->dim;
  double mod_k_effective = this->k_effective_/this->fission_test_factor_;

  std::map<std::pair<int, int>, formulation::Vector> expected_results;

  formulation::Vector expected_result_g0_a0(2), expected_result_g0_a1(2),
      expected_result_g1_a0(2), expected_result_g1_a1(2);
  expected_result_g0_a0[0] = (72.1875 + 72.1875*d)/mod_k_effective;
  expected_result_g0_a0[1] = 134.0625 * (1 + d)/mod_k_effective;
  expected_result_g0_a1[0] = (72.1875 + 144.375*d)/mod_k_effective;
  expected_result_g0_a1[1] = (134.0625 + 268.125*d)/mod_k_effective;
  expected_result_g1_a0[0] = (164.0625 + 82.03125*d)/mod_k_effective;
  expected_result_g1_a0[1] = (304.6875 + 152.34375*d)/mod_k_effective;
  expected_result_g1_a1[0] = (164.0625 * (d + 1.0))/mod_k_effective;
  expected_result_g1_a1[1] = (304.6875 * (d + 1.0))/mod_k_effective;
  expected_results.insert_or_assign({0,0}, expected_result_g0_a0);
  expected_results.insert_or_assign({1,0}, expected_result_g1_a0);
  expected_results.insert_or_assign({0,1}, expected_result_g0_a1);
  expected_results.insert_or_assign({1,1}, expected_result_g1_a1);

  for (int group = 0; group < 2; ++group) {
    for (int angle = 0; angle < 2; ++angle) {
      auto angle_it = this->quadrature_set_.begin();
      if (angle == 1)
        ++angle_it;
      auto angle_ptr = *angle_it;
      EXPECT_CALL(*this->mock_finite_element_ptr_, SetCell(this->cell_ptr_));
      EXPECT_CALL(*this->mock_finite_element_ptr_, Jacobian(_))
          .Times(2)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_finite_element_ptr_, ShapeValue(_,_))
          .Times(4)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_quadrature_set_ptr_, GetQuadraturePointIndex(_))
          .WillOnce(DoDefault());

      int out_group = !group;
      std::array<int, 3> out_index{out_group, 0, 0}, in_index{group, 0, 0};
      auto out_group_moments_map = this->out_group_moments_;
      out_group_moments_map.at(in_index) = 0;

      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(this->out_group_moments_.at(out_index)))
          .Times(1)
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(this->out_group_moments_.at(in_index)))
          .Times(0);


      system::moments::MomentVector& in_group_moment = this->group_0_moment_;
      if (group == 1)
        in_group_moment = this->group_1_moment_;

      EXPECT_CALL(*this->mock_finite_element_ptr_,
                  ValueAtQuadrature(in_group_moment))
          .Times(1)
          .WillRepeatedly(DoDefault());

      cell_vector = 0;

      EXPECT_NO_THROW({
        test_saaf.FillCellFissionSourceTerm(cell_vector, this->cell_ptr_,
                                            angle_ptr,
                                            system::EnergyGroup(group),
                                            this->k_effective_,
                                            in_group_moment,
                                            out_group_moments_map);
                      });

      std::pair<int, int> result_index{group, angle};
      cell_vector -= expected_results.at(result_index);
      EXPECT_TRUE(cell_vector.l2_norm() < 1e-10)
                << "Failed: group: " << group << " angle: " << angle;
    }
  }
}

} // namespace