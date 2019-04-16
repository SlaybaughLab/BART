#include "formulation/cfem_diffusion_stamper.h"

#include <memory>

#include "domain/tests/definition_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class CFEMDiffusionStamperTest : public ::testing::Test {
 protected:
  std::unique_ptr<domain::DefinitionMock<2>> mock_domain_ptr;
  void SetUp() override;
};

void CFEMDiffusionStamperTest::SetUp() {
  mock_domain_ptr = std::make_unique<domain::DefinitionMock<2>>();
}

TEST_F(CFEMDiffusionStamperTest, Constructor) {
  formulation::CFEM_DiffusionStamper<2> test_stamper(std::move(mock_domain_ptr));
  EXPECT_EQ(mock_domain_ptr, nullptr);
}

} // namespace