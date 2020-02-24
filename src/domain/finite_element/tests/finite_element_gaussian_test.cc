#include "domain/finite_element/finite_element_gaussian.h"

#include <vector>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <gtest/gtest.h>

#include "problem/parameter_types.h"
#include "finite_element_test.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

/* Tests for the FiniteElementGaussian class.
 *
 * As this class mostly just instantiates other classes and forwards to them
 * it's mostly interested in verifying dynamic casts to the classes we expect
 * don't return nullptr.
 */
template <typename DimensionWrapper>
class DomainFiniteElementGaussianTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  std::array<int, 4> dofs_per_cell_{1, 3, 9, 27};
  std::array<problem::DiscretizationType, 2> discretization_types {
    problem::DiscretizationType::kContinuousFEM,
    problem::DiscretizationType::kDiscontinuousFEM
  };
};

TYPED_TEST_CASE(DomainFiniteElementGaussianTest, bart::testing::AllDimensions);

// Verify that the continuous case instantiates all the correct objects
TYPED_TEST(DomainFiniteElementGaussianTest, ConstructorObjects) {
  /* Spatial dimension of the problem, this will be used to determine many of
   * the expected values for tests */
  constexpr int dim = this->dim;
  for (const auto discretization_type : this->discretization_types) {

    domain::finite_element::FiniteElementGaussian<dim> test_fe{discretization_type, 2};

    // Verify correct objects were instantiated
    ASSERT_NE(nullptr, dynamic_cast<dealii::FEValues<dim> *>(test_fe.values()));
    ASSERT_NE(nullptr,
              dynamic_cast<dealii::FEFaceValues<dim> *>(test_fe.face_values()));
    ASSERT_NE(nullptr,
              dynamic_cast<dealii::QGauss<dim> *>(test_fe.cell_quadrature()));
    ASSERT_NE(nullptr,
              dynamic_cast<dealii::QGauss<dim - 1> *>(test_fe.face_quadrature()));

    // Neighbor values are only for discontinuous so this should be a null object
    if (discretization_type == problem::DiscretizationType::kDiscontinuousFEM) {
      ASSERT_NE(nullptr,
                dynamic_cast<dealii::FE_DGQ<dim>*>(test_fe.finite_element()));
      ASSERT_NE(nullptr,
                dynamic_cast<dealii::FEFaceValues<dim> *>(
                    test_fe.neighbor_face_values()));
    } else {
      ASSERT_NE(nullptr,
                dynamic_cast<dealii::FE_Q<dim> *>(test_fe.finite_element()));
      ASSERT_EQ(nullptr,
                dynamic_cast<dealii::FEFaceValues<dim> *>(
                    test_fe.neighbor_face_values()));
    }
  }
}

// Verify that the instantiated objects were initialized correctly
TYPED_TEST(DomainFiniteElementGaussianTest, ConstructorValues) {
  /* Spatial dimension of the problem, this will be used to determine many of
   * the expected values for tests */
  constexpr int dim = this->dim;
  domain::finite_element::FiniteElementGaussian<dim> test_fe{
      problem::DiscretizationType::kContinuousFEM, 2};

  // Verify correct values returned
  ASSERT_EQ(test_fe.polynomial_degree(), 2);
  ASSERT_EQ(test_fe.dofs_per_cell(), this->dofs_per_cell_.at(dim));
  ASSERT_EQ(test_fe.n_cell_quad_pts(), this->dofs_per_cell_.at(dim));
  ASSERT_EQ(test_fe.n_face_quad_pts(), this->dofs_per_cell_.at(dim - 1));
}

TYPED_TEST(DomainFiniteElementGaussianTest, ConstructorNone) {
  constexpr int dim = this->dim;
  ASSERT_ANY_THROW({
    domain::finite_element::FiniteElementGaussian<dim> test_fe(problem::DiscretizationType::kNone, 2);
                   });
}

TYPED_TEST(DomainFiniteElementGaussianTest, ValueTest) {
  /* To use the ShapeValue, ShapeGraident, and Jacobian functions, the various
   * value objects inside our FiniteElement object need to be associated with
   * a DOF handler, which it does not own. So we will instantiate one to
   * associate it with and check that the functions are forwarding the right
   * values */
  constexpr int dim = this->dim;
  bart::domain::finite_element::FiniteElementGaussian<dim> test_fe{
      problem::DiscretizationType::kContinuousFEM, 2};

  // Triangulation and DOF handler to link to our values
  dealii::Triangulation<dim> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler dof_handler(triangulation);
  dof_handler.distribute_dofs(*test_fe.finite_element());

  test_fe.SetCell(dof_handler.begin_active());

  int cell_dofs = test_fe.dofs_per_cell();
  int cell_quad_points = test_fe.n_cell_quad_pts();

  for (int i = 0; i < cell_dofs; ++i) {
    for (int q = 0; q < cell_quad_points; ++q) {
      EXPECT_DOUBLE_EQ(test_fe.values()->shape_value(i, q),
                       test_fe.ShapeValue(i, q));
      EXPECT_TRUE(test_fe.values()->shape_grad(i, q)
                      == test_fe.ShapeGradient(i, q));
      EXPECT_DOUBLE_EQ(test_fe.values()->JxW(q),
                       test_fe.Jacobian(q));
    }
  }

  test_fe.SetFace(dof_handler.begin_active(), domain::FaceIndex(0));

  int face_quad_points = test_fe.n_face_quad_pts();
  for (int i = 0; i < cell_dofs; ++i) {
    for (int q = 0; q < face_quad_points; ++q) {
      EXPECT_DOUBLE_EQ(test_fe.face_values()->shape_value(i, q),
                       test_fe.FaceShapeValue(i, q));
      EXPECT_DOUBLE_EQ(test_fe.face_values()->JxW(q),
                       test_fe.FaceJacobian(q));
    }
  }

}

TYPED_TEST(DomainFiniteElementGaussianTest, FaceNormalTest) {
  constexpr int dim = this->dim;
  bart::domain::finite_element::FiniteElementGaussian<dim> test_fe{
      problem::DiscretizationType::kContinuousFEM, 2};

  // Triangulation and DOF handler to link to our values
  dealii::Triangulation<dim> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler dof_handler(triangulation);
  dof_handler.distribute_dofs(*test_fe.finite_element());

  test_fe.SetFace(dof_handler.begin_active(), domain::FaceIndex(0));

  EXPECT_EQ(test_fe.face_values()->normal_vector(0),
            test_fe.FaceNormal());
}

// BASE CLASS TESTS ============================================================
template <typename DimensionWrapper>
class DomainFiniteElementGaussianBaseMethodsTest :
    public domain::finite_element::testing::FiniteElementBaseClassTest<DimensionWrapper::value> {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  void SetUp() override {
    domain::finite_element::testing::FiniteElementBaseClassTest<dim>::SetUp();
  }
};

TYPED_TEST_CASE(DomainFiniteElementGaussianBaseMethodsTest,
                bart::testing::AllDimensions);

TYPED_TEST(DomainFiniteElementGaussianBaseMethodsTest, BaseSetCell) {
  bart::domain::finite_element::FiniteElementGaussian<this->dim> test_fe{
    problem::DiscretizationType::kDiscontinuousFEM, 2};
  this->TestSetCell(&test_fe);
}

TYPED_TEST(DomainFiniteElementGaussianBaseMethodsTest, BaseSetCellAndFace) {
  bart::domain::finite_element::FiniteElementGaussian<this->dim> test_fe{
      problem::DiscretizationType::kDiscontinuousFEM, 2};
  this->TestSetCellAndFace(&test_fe);
}

TYPED_TEST(DomainFiniteElementGaussianBaseMethodsTest, BaseValueAtQuadrature) {
  bart::domain::finite_element::FiniteElementGaussian<this->dim> test_fe{
      problem::DiscretizationType::kDiscontinuousFEM, 2};
  this->TestValueAtQuadrature(&test_fe);
}

} // namespace