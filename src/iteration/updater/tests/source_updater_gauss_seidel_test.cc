#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "test_helpers/gmock_wrapper.h"
#include "data/system/system.h"
#include "data/system/system_types.h"
#include "data/system/tests/right_hand_side_mock.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/cfem_stamper_i.h"


namespace  {

using namespace bart;
using data::system::MPIVector;

using ::testing::An;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::_;

class IterationSourceUpdaterGaussSeidelTest : public ::testing::Test {
 protected:
  using CFEMSourceUpdater = iteration::updater::SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;

  std::unique_ptr<formulation::CFEM_StamperMock> mock_stamper_ptr_;
  std::unique_ptr<data::system::RightHandSideMock> mock_rhs_ptr_;
  data::system::System test_system_;

  std::shared_ptr<MPIVector> source_vector_ptr_;

  void SetUp() override;
};

void IterationSourceUpdaterGaussSeidelTest::SetUp() {
  mock_stamper_ptr_ = std::make_unique<formulation::CFEM_StamperMock>();
  mock_rhs_ptr_ = std::make_unique<NiceMock<data::system::RightHandSideMock>>();

  ON_CALL(*mock_rhs_ptr_, GetVariablePtr(An<data::system::Index>(),_))
      .WillByDefault(Return(source_vector_ptr_));
  ON_CALL(*mock_rhs_ptr_, GetVariablePtr(An<data::system::GroupNumber>(),_))
      .WillByDefault(Return(source_vector_ptr_));

  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
}

TEST_F(IterationSourceUpdaterGaussSeidelTest, Constructor) {
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
  EXPECT_EQ(mock_stamper_ptr_, nullptr);
}

} // namespace