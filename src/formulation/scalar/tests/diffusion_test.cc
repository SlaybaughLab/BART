#include "formulation/scalar/cfem_diffusion.h"

#include <memory>

#include "data/cross_sections.h"
#include "domain/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"



namespace {

using ::testing::DoDefault;
using ::testing::NiceMock;
using ::testing::Return;
using namespace bart;

class FormulationCFEMDiffusionTest : public ::testing::Test {
 protected:
  std::shared_ptr<domain::FiniteElementMock<2>> fe_mock_ptr;
  std::shared_ptr<data::CrossSections> cross_sections_ptr;
  void SetUp() override;
};

void FormulationCFEMDiffusionTest::SetUp() {
  // Make mock objects. Cross-sections is a struct that cannot be mocked, but
  // we can mock the material object it is based on.
  NiceMock<btest::MockMaterial> mock_material;
  fe_mock_ptr = std::make_shared<NiceMock<domain::FiniteElementMock<2>>>();
  cross_sections_ptr = std::make_shared<data::CrossSections>(mock_material);

  ON_CALL(*fe_mock_ptr, dofs_per_cell())
      .WillByDefault(Return(4));
  ON_CALL(*fe_mock_ptr, n_cell_quad_pts())
      .WillByDefault(Return(2));
  ON_CALL(*fe_mock_ptr, n_face_quad_pts())
      .WillByDefault(Return(2));
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

}

} // namespace