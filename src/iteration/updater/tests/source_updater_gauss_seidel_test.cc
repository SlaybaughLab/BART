#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"
#include "data/moment_types.h"
#include "data/system.h"
#include "data/system/system_types.h"
#include "data/system/tests/right_hand_side_mock.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/cfem_stamper_i.h"


namespace  {

using namespace bart;
using data::system::MPIVector;

using ::testing::An;
using ::testing::Ref;
using ::testing::DoDefault;
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
  using VariableTerms = data::system::RightHandSideI::VariableTerms;

  // Required objects
  // Test system
  data::System test_system_;
  // Vector returned by the mock RightHandSide object when the variable
  // right hand side vector is requested. This is what should be stamped.
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

  data::system::MomentsMap current_iteration, previous_iteration;

  for (data::system::GroupNumber group = 0; group < 5; ++group) {
    for (data::system::HarmonicL l = 0; l < 2; ++l) {
      for (data::system::HarmonicM m = -l; m <= l; ++m) {
        data::system::MomentVector current_moment, previous_moment;
        current_iteration[{group, l, m}] = current_moment;
        previous_iteration[{group, l, m}] = previous_moment;
      }
    }
  }
}

// Verifies that the Updater takes ownership of the stamper.
TEST_F(IterationSourceUpdaterGaussSeidelTest, Constructor) {
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);

  EXPECT_EQ(mock_stamper_ptr_, nullptr);
}

// Verifies operation of the UpdateScatteringSource function
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateScatteringSourceTest) {
  data::system::GroupNumber group = btest::RandomDouble(0, 6);
  data::system::AngleIndex angle = btest::RandomDouble(0, 10);

  EXPECT_CALL(*mock_rhs_ptr_, GetVariablePtr(group, VariableTerms::kScatteringSource))
      .WillOnce(DoDefault());
  EXPECT_CALL(*mock_stamper_ptr_,
      StampScatteringSource(Ref(*source_vector_ptr_), // Vector to stamp from mock RHS
                            group,                    // Group specified by the test
                            Ref(test_system_.current_iteration[{group, 0, 0}]), // Current scalar flux for in-group
                            Ref(test_system_.current_iteration)));              // Current moments for out-group

  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
  test_updater.UpdateScatteringSource(test_system_, group, angle);


}

} // namespace