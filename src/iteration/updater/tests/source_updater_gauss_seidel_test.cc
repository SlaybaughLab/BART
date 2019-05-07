#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "test_helpers/gmock_wrapper.h"
#include "data/moment_types.h"
#include "data/system/system.h"
#include "data/system/system_types.h"
#include "data/system/tests/right_hand_side_mock.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/cfem_stamper_i.h"


namespace  {

using namespace bart;
using data::system::MPIVector;

using ::testing::An;
using ::testing::ByRef;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::_;

/* This test class verifies the operation of the SourceUpdaterGaussSeidel class.
 * It uses the template of it that uses a CFEMStamperI, a mock version of
 * that object will be used. The required System object will be explicitly
 * created, because it is just a struct, but some of the stored objects
 * (the right hand side) will also use mocks.
 */
class IterationSourceUpdaterGaussSeidelTest : public ::testing::Test {
 protected:
  using CFEMSourceUpdater = iteration::updater::SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;

  // Required objects
  data::system::System test_system_;
  std::shared_ptr<MPIVector> source_vector_ptr_;

  // Required mocks
  std::unique_ptr<formulation::CFEM_StamperMock> mock_stamper_ptr_;
  std::unique_ptr<data::system::RightHandSideMock> mock_rhs_ptr_;

  void SetUp() override;
};

/* Sets up the tests. This function first creates the mock objects to be used
 * by the testing, then establishes any default behaviors for NiceMocks. The
 * right hand side vector can be handed off to the System object, but the
 * stamper needs to be given to the test updater _in_ the tests because, as a
 * required dependency, it is passed in the constructor.
 */
void IterationSourceUpdaterGaussSeidelTest::SetUp() {
  mock_stamper_ptr_ = std::make_unique<formulation::CFEM_StamperMock>();
  mock_rhs_ptr_ = std::make_unique<NiceMock<data::system::RightHandSideMock>>();

  ON_CALL(*mock_rhs_ptr_, GetVariablePtr(An<data::system::Index>(),_))
      .WillByDefault(Return(source_vector_ptr_));
  ON_CALL(*mock_rhs_ptr_, GetVariablePtr(An<data::system::GroupNumber>(),_))
      .WillByDefault(Return(source_vector_ptr_));

  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
}

// Verifies that the Updater takes ownership of the stamper.
TEST_F(IterationSourceUpdaterGaussSeidelTest, Constructor) {
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
  EXPECT_EQ(mock_stamper_ptr_, nullptr);
}

// Verifies operation of the UpdateScatteringSource function
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateScatteringSourceTest) {
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
}

} // namespace