#include "formulation/scalar/cfem_diffusion.h"

#include <memory>

#include "data/cross_sections.h"
#include "domain/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"



namespace {

using ::testing::NiceMock;
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
  fe_mock_ptr = std::make_shared<domain::FiniteElementMock<2>>();
  cross_sections_ptr = std::make_shared<data::CrossSections>(mock_material);
}

TEST_F(FormulationCFEMDiffusionTest, ConstructorTest) {
  formulation::scalar::CFEM_Diffusion<2> test_diffusion(fe_mock_ptr,
                                                        cross_sections_ptr);
}

} // namespace