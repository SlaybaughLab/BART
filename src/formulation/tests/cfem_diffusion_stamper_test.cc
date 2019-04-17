#include "formulation/cfem_diffusion_stamper.h"

#include <memory>

#include <deal.II/lac/full_matrix.h>

#include "domain/tests/definition_mock.h"
#include "formulation/scalar/tests/cfem_diffusion_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using ::testing::DoDefault;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::Unused;
using ::testing::_;

using namespace bart;
using Cell = domain::DefinitionI<2>::Cell;
using InitToken = formulation::scalar::CFEM_DiffusionI<2>::InitializationToken;
using Matrix = dealii::FullMatrix<double>;

class CFEMDiffusionStamperTest : public ::testing::Test {
 protected:
  std::unique_ptr<NiceMock<domain::DefinitionMock<2>>> mock_definition_ptr;
  std::unique_ptr<NiceMock<formulation::scalar::CFEM_DiffusionMock<2>>> mock_diffusion_ptr;
  void SetUp() override;
  InitToken init_token_;
};

// TODO(Josh) Move this to a header where other stamper tests can use it
void FillMatrixWithOnes(Matrix& to_fill, Unused, Unused, Unused, Unused) {
  for (int i = 0; i < to_fill.n_rows(); ++i) {
    for (int j = 0; j < to_fill.n_cols(); ++j) {
      to_fill(i,j) = 1;
    }
  }
}

void CFEMDiffusionStamperTest::SetUp() {
  mock_definition_ptr = std::make_unique<NiceMock<domain::DefinitionMock<2>>>();
  mock_diffusion_ptr =
      std::make_unique<NiceMock<formulation::scalar::CFEM_DiffusionMock<2>>>();

  ON_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillByDefault(Return(init_token_));
  ON_CALL(*mock_diffusion_ptr, FillCellStreamingTerm(_, _, _, _, _))
      .WillByDefault(Invoke(FillMatrixWithOnes));
}

TEST_F(CFEMDiffusionStamperTest, Constructor) {
  Cell test_cell;
  std::vector<Cell> cells{test_cell};

  EXPECT_CALL(*mock_definition_ptr, Cells())
      .WillOnce(Return(cells));
  EXPECT_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  EXPECT_EQ(mock_diffusion_ptr, nullptr);
  EXPECT_EQ(mock_definition_ptr, nullptr);
}

TEST_F(CFEMDiffusionStamperTest, FillCellStreamingTest) {

}



} // namespace