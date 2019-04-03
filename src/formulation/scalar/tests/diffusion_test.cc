#include "formulation/scalar/cfem_diffusion.h"

#include <array>
#include <memory>
#include <cstdlib>

#include "data/cross_sections.h"
#include "domain/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"


namespace {

using ::testing::AssertionResult;
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;
using ::testing::Le;
using ::testing::Ge;
using ::testing::DoDefault;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::_;
using namespace bart;

class FormulationCFEMDiffusionTest : public ::testing::Test {
 protected:
  using Matrix = dealii::FullMatrix<double>;
  std::shared_ptr<domain::FiniteElementMock<2>> fe_mock_ptr;
  std::shared_ptr<data::CrossSections> cross_sections_ptr;

  int *a, *b;

  void SetUp() override;
};

void FormulationCFEMDiffusionTest::SetUp() {
  // Make mock objects. Cross-sections is a struct that cannot be mocked, but
  // we can mock the material object it is based on.
  NiceMock<btest::MockMaterial> mock_material;
  fe_mock_ptr = std::make_shared<NiceMock<domain::FiniteElementMock<2>>>();
  cross_sections_ptr = std::make_shared<data::CrossSections>(mock_material);

  ON_CALL(*fe_mock_ptr, dofs_per_cell())
      .WillByDefault(Return(2));
  ON_CALL(*fe_mock_ptr, n_cell_quad_pts())
      .WillByDefault(Return(2));
  ON_CALL(*fe_mock_ptr, n_face_quad_pts())
      .WillByDefault(Return(2));
}

AssertionResult CompareMatrices(const dealii::FullMatrix<double>& expected,
                                const dealii::FullMatrix<double>& result,
                                const double tol = 1e-6) {
  unsigned int rows = expected.m();
  unsigned int cols = expected.n();

  if (result.m() != rows)
    return AssertionFailure() << "Result has wrong number of rows: "
                              << result.m() << ", expected" << rows;
  if (result.n() != cols)
    return AssertionFailure() << "Result has wrong number of columns: "
                              << result.n() << ", expected" << cols;

  for (unsigned int i = 0; i < rows; ++i) {
    for (unsigned int j = 0; j < cols; ++j) {
      if (abs(result(i, j) - expected(i, j)) > tol) {
        return AssertionFailure() << "Entry (" << i << ", " << j <<
                                  ") has value: " << result(i, j) <<
                                  ", expected: " << expected(i, j);
      }
    }
  }
  return AssertionSuccess();
}


TEST_F(FormulationCFEMDiffusionTest, ConstructorTest) {
  EXPECT_CALL(*fe_mock_ptr, dofs_per_cell())
      .WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, n_cell_quad_pts())
      .WillOnce(DoDefault());
  EXPECT_CALL(*fe_mock_ptr, n_face_quad_pts())
      .WillOnce(DoDefault());

  formulation::scalar::CFEM_Diffusion<2> test_diffusion(fe_mock_ptr,
                                                        cross_sections_ptr);

  EXPECT_EQ(fe_mock_ptr.use_count(), 2);
  EXPECT_EQ(cross_sections_ptr.use_count(), 2);
}

TEST_F(FormulationCFEMDiffusionTest, PrecalculateShapeTest) {
  // Shape call, by default, returns (quadrature point index + degree of freedom)

  for (int q = 0; q < 2; ++q) {
    for (int i = 0; i < 2; ++i) {
      EXPECT_CALL(*fe_mock_ptr, ShapeValue(i, q))
      .Times(4)
      .WillRepeatedly(Return(i + q));
    }
  }

  std::array<double, 4> shape_matrix_q_0_values = {0, 0,
                                                   0, 1};
  std::array<double, 4> shape_matrix_q_1_values = {1, 2,
                                                   2, 4};

  dealii::FullMatrix<double> shape_matrix_q_0{2,2, shape_matrix_q_0_values.begin()};
  dealii::FullMatrix<double> shape_matrix_q_1{2,2, shape_matrix_q_1_values.begin()};

  formulation::scalar::CFEM_Diffusion<2> test_diffusion(fe_mock_ptr,
                                                        cross_sections_ptr);
  dealii::DoFHandler<2>::active_cell_iterator it;

  test_diffusion.Precalculate(it);
  auto shape_squared = test_diffusion.GetShapeSquared();

  EXPECT_TRUE(CompareMatrices(shape_matrix_q_0, shape_squared.at(0)));
  EXPECT_TRUE(CompareMatrices(shape_matrix_q_1, shape_squared.at(1)));
}



} // namespace