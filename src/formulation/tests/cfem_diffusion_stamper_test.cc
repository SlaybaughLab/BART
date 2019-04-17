#include "formulation/cfem_diffusion_stamper.h"

#include <memory>

#include "domain/tests/definition_mock.h"
#include "formulation/scalar/tests/cfem_diffusion_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using ::testing::Return;
using ::testing::_;
using namespace bart;
using Cell = domain::DefinitionI<2>::Cell;
using InitToken = formulation::scalar::CFEM_DiffusionI<2>::InitializationToken;

class CFEMDiffusionStamperTest : public ::testing::Test {
 protected:
  std::unique_ptr<domain::DefinitionMock<2>> mock_definition_ptr;
  std::unique_ptr<formulation::scalar::CFEM_DiffusionMock<2>> mock_diffusion_ptr;
  void SetUp() override;
};

void CFEMDiffusionStamperTest::SetUp() {
  mock_definition_ptr = std::make_unique<domain::DefinitionMock<2>>();
  mock_diffusion_ptr =
      std::make_unique<formulation::scalar::CFEM_DiffusionMock<2>>();
}

TEST_F(CFEMDiffusionStamperTest, Constructor) {
  InitToken init_token;
  Cell test_cell;
  std::vector<Cell> cells{test_cell};

  EXPECT_CALL(*mock_definition_ptr, Cells())
      .WillOnce(Return(cells));
  EXPECT_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillOnce(Return(init_token));

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  EXPECT_EQ(mock_diffusion_ptr, nullptr);
  EXPECT_EQ(mock_definition_ptr, nullptr);
}



} // namespace