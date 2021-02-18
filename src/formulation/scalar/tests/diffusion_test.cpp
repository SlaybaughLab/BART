#include "formulation/scalar/diffusion.hpp"

#include <array>
#include <memory>
#include <cstdlib>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "data/cross_sections/material_cross_sections.hpp"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "data/material/tests/material_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"


namespace {

using ::testing::Le;
using ::testing::Ge;
using ::testing::DoDefault;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::_;
using ::testing::Ref;
using ::testing::AtLeast;
using namespace bart;

namespace test_helpers = bart::test_helpers;
using test_helpers::AreEqual;

class FormulationCFEMDiffusionTest : public ::testing::Test {
 protected:
  FormulationCFEMDiffusionTest()
      : dof_handler_(triangulation_),
        fe_(1) {};
  using Matrix = dealii::FullMatrix<double>;
  std::shared_ptr<domain::finite_element::FiniteElementMock<2>> fe_mock_ptr;
  std::shared_ptr<data::cross_sections::MaterialCrossSections> cross_sections_ptr;

  dealii::DoFHandler<2>::active_cell_iterator cell_ptr_;
  dealii::Triangulation<2> triangulation_;
  dealii::DoFHandler<2> dof_handler_;
  dealii::FE_Q<2> fe_;

  const int fissile_material_id_ = 0, non_fissile_material_id_ = 1;

  void SetUp() override;
  void SetUpDealii();
};

void FormulationCFEMDiffusionTest::SetUp() {
  SetUpDealii();
  // Make mock objects. Cross-sections is a struct that cannot be mocked, but
  // we can mock the material object it is based on.
  NiceMock<data::material::MaterialMock> mock_material;
  fe_mock_ptr = std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<2>>>();

  ON_CALL(*fe_mock_ptr, dofs_per_cell()).WillByDefault(Return(2));
  ON_CALL(*fe_mock_ptr, n_cell_quad_pts()).WillByDefault(Return(2));
  ON_CALL(*fe_mock_ptr, n_face_quad_pts()).WillByDefault(Return(2));

  for (int q = 0; q < 2; ++q) {
    ON_CALL(*fe_mock_ptr, Jacobian(q)).WillByDefault(Return((1 + q)*3));
    ON_CALL(*fe_mock_ptr, FaceJacobian(q)).WillByDefault(Return((1 + q)*3));
    for (int i = 0; i < 2; ++i) {
      ON_CALL(*fe_mock_ptr, ShapeValue(i, q)).WillByDefault(Return(i + q));
      ON_CALL(*fe_mock_ptr, FaceShapeValue(i,q)).WillByDefault(Return(i + q));

      dealii::Tensor<1, 2> gradient_tensor;
      gradient_tensor[0] = i;
      gradient_tensor[1] = q;

      ON_CALL(*fe_mock_ptr, ShapeGradient(i, q)).WillByDefault(Return(gradient_tensor));
    }
  }

  // Cross-section data for fake two-group
  std::array<double, 4> sigma_s_values{0.25, 0.5, 0.75, 1.0};
  dealii::FullMatrix<double> sigma_s_matrix{2,2, sigma_s_values.begin()};
  std::unordered_map<int, dealii::FullMatrix<double>> sigma_s{{this->fissile_material_id_, sigma_s_matrix}};

  std::unordered_map<int, std::vector<double>> sigma_t{{this->fissile_material_id_, {1.0, 2.0}}}; // i + 1

  std::unordered_map<int, std::vector<double>> q{{this->non_fissile_material_id_, {1.0, 2.0}}};

  std::unordered_map<int, bool> fissile_id{{this->fissile_material_id_, true}, 
                                           {this->non_fissile_material_id_, false}};

  ON_CALL(mock_material, GetSigT()).WillByDefault(Return(sigma_t));
  ON_CALL(mock_material, GetSigS()).WillByDefault(Return(sigma_s));
  ON_CALL(mock_material, GetDiffusionCoef()).WillByDefault(Return(sigma_t));
  ON_CALL(mock_material, GetQ()).WillByDefault(Return(q));
  ON_CALL(mock_material, GetFissileIDMap()).WillByDefault(Return(fissile_id));
  ON_CALL(mock_material, GetChiNuSigF()).WillByDefault(Return(sigma_s));

  cross_sections_ptr = std::make_shared<data::cross_sections::MaterialCrossSections>(mock_material);
}

// Set up a simple deal.ii problem so that the cell points to something, this
// is for material ID
void FormulationCFEMDiffusionTest::SetUpDealii() {
  dealii::GridGenerator::hyper_cube(triangulation_, 0, 1);
  dof_handler_.distribute_dofs(fe_);
  for (auto cell = dof_handler_.begin_active(); cell != dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned()) {
      cell_ptr_ = cell;
      cell_ptr_->set_material_id(this->fissile_material_id_);
    }
  }
}

TEST_F(FormulationCFEMDiffusionTest, ConstructorTest) {
  EXPECT_CALL(*fe_mock_ptr, dofs_per_cell()).WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, n_cell_quad_pts()).WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, n_face_quad_pts()).WillOnce(DoDefault());

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr,cross_sections_ptr);

  EXPECT_EQ(fe_mock_ptr.use_count(), 3);
  EXPECT_EQ(cross_sections_ptr.use_count(), 2);
  EXPECT_FALSE(test_diffusion.is_initialized());
}

TEST_F(FormulationCFEMDiffusionTest, PrecalculateTest) {
  // Shape call, by default, returns (quadrature point index + degree of freedom)
  // Gradient call will return a tensor of [i, q]
  // Expected results, we need to make arrays to put them into dealii matrices
  std::array<double, 4> shape_matrix_q_0_values = {0, 0,
                                                   0, 1};
  std::array<double, 4> shape_matrix_q_1_values = {1, 2,
                                                   2, 4};
  std::array<double, 4> gradient_matrix_q_0_values = {0, 0,
                                                      0, 1};
  std::array<double, 4> gradient_matrix_q_1_values = {1, 1,
                                                      1, 2};

  dealii::FullMatrix<double> shape_matrix_q_0{2,2,shape_matrix_q_0_values.begin()};
  dealii::FullMatrix<double> shape_matrix_q_1{2,2,shape_matrix_q_1_values.begin()};
  dealii::FullMatrix<double> gradient_matrix_q_0{2,2,gradient_matrix_q_0_values.begin()};
  dealii::FullMatrix<double> gradient_matrix_q_1{2,2,gradient_matrix_q_1_values.begin()};

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr,cross_sections_ptr);

  // Set call expectations
  EXPECT_CALL(*fe_mock_ptr, SetCell(_)).Times(1);
  // These are called at least 4 times which is the minimum required for 2 dofs and 2 quadrature points
  EXPECT_CALL(*fe_mock_ptr, ShapeValue(_,_)).Times(::testing::AtLeast(4)).WillRepeatedly(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, ShapeGradient(_,_)).Times(::testing::AtLeast(4)).WillRepeatedly(DoDefault());

  test_diffusion.Precalculate(cell_ptr_);
  auto shape_squared = test_diffusion.GetShapeSquared();
  auto gradient_squared = test_diffusion.GetGradientSquared();

  EXPECT_TRUE(AreEqual(shape_matrix_q_0, shape_squared.at(0)));
  EXPECT_TRUE(AreEqual(shape_matrix_q_1, shape_squared.at(1)));
  EXPECT_TRUE(AreEqual(gradient_matrix_q_0, gradient_squared.at(0)));
  EXPECT_TRUE(AreEqual(gradient_matrix_q_1, gradient_squared.at(1)));
  EXPECT_TRUE(test_diffusion.is_initialized());
}

TEST_F(FormulationCFEMDiffusionTest, FillCellConstantTermTest) {
  dealii::Vector<double> test_vector(2);
  dealii::Vector<double> expected_values{ 54, 129 };
  dealii::Vector<double> constant_vector_at_dofs(2);
  std::vector<double> constant_vector_at_quadrature{ 7, 9 };

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);
  EXPECT_CALL(*fe_mock_ptr, dofs_per_cell()).WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, n_cell_quad_pts()).WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, SetCell(cell_ptr_)).Times(1);
  EXPECT_CALL(*fe_mock_ptr, ValueAtQuadrature(Ref(constant_vector_at_dofs)))
      .WillOnce(Return(constant_vector_at_quadrature));

  for (int q = 0; q < 2; ++q) {
    EXPECT_CALL(*fe_mock_ptr, Jacobian(q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
    for (int i = 0; i < 2; ++i) {
      EXPECT_CALL(*fe_mock_ptr, ShapeValue(i,q)).Times(AtLeast(1)).WillRepeatedly(DoDefault());
    }
  }

  test_diffusion.FillCellConstantTerm(test_vector, cell_ptr_, constant_vector_at_dofs);
  EXPECT_TRUE(AreEqual(expected_values, test_vector));
}

TEST_F(FormulationCFEMDiffusionTest, FillCellStreamingTermTest) {
  dealii::FullMatrix<double> test_matrix(2,2);

  std::array<double, 4> expected_values{6, 6,
                                        6, 15};
  dealii::FullMatrix<double> expected_matrix(2, 2, expected_values.begin());

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  EXPECT_CALL(*fe_mock_ptr, Jacobian(_)).Times(2).WillRepeatedly(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, SetCell(cell_ptr_)).Times(2);

  EXPECT_ANY_THROW({
    test_diffusion.FillCellStreamingTerm(test_matrix, cell_ptr_, 0);
  });
  EXPECT_NO_THROW({
    test_diffusion.Precalculate(cell_ptr_);
    test_diffusion.FillCellStreamingTerm(test_matrix, cell_ptr_, 0);
  });

  EXPECT_TRUE(AreEqual(expected_matrix, test_matrix));
}

TEST_F(FormulationCFEMDiffusionTest, FillCellCollisionTermTest) {
  dealii::FullMatrix<double> test_matrix(2,2);

  std::array<double, 4> expected_values{4.5, 9.0,
                                        9.0, 20.25};
  dealii::FullMatrix<double> expected_matrix(2, 2, expected_values.begin());

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  EXPECT_CALL(*fe_mock_ptr, Jacobian(_)).Times(2).WillRepeatedly(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, SetCell(cell_ptr_)).Times(2);

  EXPECT_ANY_THROW({
    test_diffusion.FillCellCollisionTerm(test_matrix, cell_ptr_, 0);
                   });
  EXPECT_NO_THROW({
    test_diffusion.Precalculate(cell_ptr_);
    test_diffusion.FillCellCollisionTerm(test_matrix, cell_ptr_, 0);
                  });

  EXPECT_TRUE(AreEqual(expected_matrix, test_matrix));
}

TEST_F(FormulationCFEMDiffusionTest, FillBoundaryTermTestReflective) {
  dealii::FullMatrix<double> test_matrix(2,2);
  dealii::FullMatrix<double> expected_matrix(2,2);

  auto boundary = formulation::scalar::Diffusion<2>::BoundaryType::kReflective;

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  EXPECT_ANY_THROW({
    test_diffusion.FillBoundaryTerm(test_matrix, cell_ptr_, 0, boundary);
                   });
  EXPECT_NO_THROW({
    test_diffusion.Precalculate(cell_ptr_);
    test_diffusion.FillBoundaryTerm(test_matrix, cell_ptr_, 0, boundary);
                  });

  EXPECT_TRUE(AreEqual(expected_matrix, test_matrix));
}

TEST_F(FormulationCFEMDiffusionTest, FillBoundaryTermTestVacuum) {
  dealii::FullMatrix<double> test_matrix(2,2);

  auto boundary = formulation::scalar::Diffusion<2>::BoundaryType::kVacuum;

  std::array<double, 4> expected_values{3, 6,
                                        6, 13.5};

  dealii::FullMatrix<double> expected_matrix(2, 2, expected_values.begin());

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  EXPECT_CALL(*fe_mock_ptr, SetFace(cell_ptr_, domain::FaceIndex(0))).Times(1);
  EXPECT_CALL(*fe_mock_ptr, FaceJacobian(_)).Times(2).WillRepeatedly(DoDefault());

  test_diffusion.Precalculate(cell_ptr_);
  test_diffusion.FillBoundaryTerm(test_matrix, cell_ptr_, 0, boundary);

  EXPECT_TRUE(AreEqual(expected_matrix, test_matrix));

}

TEST_F(FormulationCFEMDiffusionTest, FillCellFixedSourceNonFissile) {
  dealii::Vector<double> expected_vector{6, 15};
  dealii::Vector<double> test_vector(2);

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr,cross_sections_ptr);

  cell_ptr_->set_material_id(this->non_fissile_material_id_);
  EXPECT_CALL(*fe_mock_ptr, SetCell(_)).Times(1);
  EXPECT_CALL(*fe_mock_ptr, Jacobian(_)).Times(2).WillRepeatedly(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, ShapeValue(_,_)).Times(4).WillRepeatedly(DoDefault());

  EXPECT_NO_THROW(test_diffusion.FillCellFixedSource(test_vector, cell_ptr_, 0));

  EXPECT_TRUE(AreEqual(expected_vector, test_vector));
}

TEST_F(FormulationCFEMDiffusionTest, FillCellFixedSourceFissile) {
  dealii::Vector<double> expected_vector{0, 0};
  dealii::Vector<double> test_vector(2);

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  EXPECT_NO_THROW(test_diffusion.FillCellFixedSource(test_vector, cell_ptr_, 0));

  EXPECT_TRUE(AreEqual(expected_vector, test_vector));
}



TEST_F(FormulationCFEMDiffusionTest, FillFissionSourceTest) {

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  dealii::Vector<double> test_vector(2);
  int group = 0;
  double k_effective = 1.05;
  // Make in-group moment
  std::vector<double> in_group_moment_values{0.5, 0.5};
  system::moments::MomentVector in_group_moment(in_group_moment_values.begin(), in_group_moment_values.end());

  /* Make out-group moments (specifically in-group values in this are different
   * This is a somewhat tortured process due to the way dealii::Vectors are
   * defined (system::moments::MomentVector is an alias) */
  std::vector<double> group_1_moment_values{1.0, 1.0};
  system::moments::MomentVector group_0_moment{0.75, 0.75};
  system::moments::MomentVector group_1_moment{1.0, 1.0};
  system::moments::MomentsMap out_group_moments;
  out_group_moments[{0,0,0}] = group_0_moment;
  out_group_moments[{1,0,0}] = group_1_moment;

  // Expected answer
  dealii::Vector<double> expected_vector{5.0, 12.5};

  EXPECT_CALL(*fe_mock_ptr, SetCell(_)).Times(1);
  EXPECT_CALL(*fe_mock_ptr, ValueAtQuadrature(group_1_moment)).WillOnce(Return(group_1_moment_values));
  EXPECT_CALL(*fe_mock_ptr, ValueAtQuadrature(in_group_moment)).WillOnce(Return(in_group_moment_values));

  test_diffusion.FillCellFissionSource(test_vector, cell_ptr_, group, k_effective, in_group_moment, out_group_moments);

  EXPECT_TRUE(AreEqual(expected_vector, test_vector));

}

TEST_F(FormulationCFEMDiffusionTest, FillScatteringSourceTest) {

  formulation::scalar::Diffusion<2> test_diffusion(fe_mock_ptr, cross_sections_ptr);

  dealii::Vector<double> test_vector(2);
  // Make in-group moment
  std::vector<double> in_group_moment_values{0.5, 0.5};
  system::moments::MomentVector in_group_moment(in_group_moment_values.begin(), in_group_moment_values.end());

  /* Make out-group moments (specifically in-group values in this are different
   * This is a somewhat tortured process due to the way dealii::Vectors are
   * defined (system::moments::MomentVector is an alias) */
  std::vector<double> group_0_moment_values{0.75, 0.75};
  std::vector<double> group_1_moment_values{1.0, 1.0};
  system::moments::MomentVector group_0_moment(group_0_moment_values.begin(), group_0_moment_values.end());
  system::moments::MomentVector group_1_moment(group_1_moment_values.begin(), group_1_moment_values.end());
  system::moments::MomentsMap out_group_moments;
  out_group_moments[{0,0,0}] = group_0_moment;
  out_group_moments[{1,0,0}] = group_1_moment;

  // Expected solutions
  std::vector<double> expected_group_0_values{3.0, 7.5};
  dealii::Vector<double> expected_group_0_vector(expected_group_0_values.begin(), expected_group_0_values.end());
  std::vector<double> expected_group_1_values{3.375, 8.4375};
  dealii::Vector<double> expected_group_1_vector(expected_group_1_values.begin(), expected_group_1_values.end());
  std::array<dealii::Vector<double>, 2> expected_vectors{ expected_group_0_vector, expected_group_1_vector };

  EXPECT_CALL(*fe_mock_ptr, SetCell(_)).Times(2);
  EXPECT_CALL(*fe_mock_ptr, ValueAtQuadrature(group_1_moment)).WillOnce(Return(group_1_moment_values));
  EXPECT_CALL(*fe_mock_ptr, ValueAtQuadrature(group_0_moment)).WillOnce(Return(group_0_moment_values));

  for (int group = 0; group < 2; ++group) {
    test_diffusion.FillCellScatteringSource(test_vector, cell_ptr_, group, out_group_moments);

    EXPECT_TRUE(AreEqual(expected_vectors.at(group), test_vector));
    test_vector = 0;
  }
}


} // namespace