#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/system.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/system_types.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/cfem_stamper_i.h"
#include "test_helpers/test_assertions.h"
#include "system/terms/tests/linear_term_mock.h"


namespace  {

using namespace bart;
using system::MPIVector;

using ::testing::An;
using ::testing::Ref;
using ::testing::DoDefault;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Return, ::testing::ReturnRef;
using ::testing::WithArgs;
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
  using VariableTerms = system::terms::VariableLinearTerms;

  // Required objects
  // Test system
  system::System test_system_;
  // Vector returned by the mock RightHandSide object when the variable
  // right hand side vector is requested. This is what should be stamped.
  std::shared_ptr<MPIVector> source_vector_ptr_;
  MPIVector expected_vector_;

  // Required mocks and supporting objects
  std::shared_ptr<formulation::CFEM_StamperMock> mock_stamper_ptr_;
  std::unique_ptr<system::terms::LinearTermMock> mock_rhs_ptr_;
  std::unique_ptr<system::moments::SphericalHarmonicMock> moments_ptr_;

  system::moments::SphericalHarmonicMock* moments_obs_ptr_;

  system::moments::MomentsMap current_iteration_moments_,
      previous_iteration_moments_;

  void SetUp() override;
};

/* Sets up the tests. This function first creates the mock objects to be used
 * by the testing, then establishes any default behaviors for NiceMocks. The
 * right hand side vector can be handed off to the System object, but the
 * stamper needs to be given to the test updater _in_ the tests because, as a
 * required dependency, it is passed in the constructor.
 */
void IterationSourceUpdaterGaussSeidelTest::SetUp() {
  mock_stamper_ptr_ = std::make_shared<NiceMock<formulation::CFEM_StamperMock>>();
  mock_rhs_ptr_ = std::make_unique<NiceMock<system::terms::LinearTermMock>>();
  moments_ptr_ = std::make_unique<system::moments::SphericalHarmonicMock>();

  /* Create and populate moment maps. The inserted MomentVectors can be empty
   * because we will check that the correct ones are passed by reference not
   * entries.
   */

  int l_max = 2;

  for (system::GroupNumber group = 0; group < 5; ++group) {
    for (system::moments::HarmonicL l = 0; l < l_max; ++l) {
      for (system::moments::HarmonicM m = -l_max; m <= l_max; ++m) {
        system::moments::MomentVector current_moment, previous_moment;
        current_iteration_moments_[{group, l, m}] = current_moment;
        previous_iteration_moments_[{group, l, m}] = previous_moment;
      }
    }
  }

  moments_obs_ptr_ = moments_ptr_.get();
  test_system_.current_moments = std::move(moments_ptr_);

  /* Initialize MPI Vectors */
  source_vector_ptr_ = std::make_shared<MPIVector>();
  auto n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  source_vector_ptr_->reinit(MPI_COMM_WORLD,
                             n_processes*5,
                             5);
  expected_vector_.reinit(*source_vector_ptr_);

  ON_CALL(*mock_rhs_ptr_, GetVariableTermPtr(An<system::Index>(),_))
      .WillByDefault(Return(source_vector_ptr_));
  ON_CALL(*mock_rhs_ptr_, GetVariableTermPtr(An<system::GroupNumber>(),_))
      .WillByDefault(Return(source_vector_ptr_));
}
// Fills an MPI vector with value
void StampMPIVector(MPIVector &to_fill, double value = 2) {
  auto [local_begin, local_end] = to_fill.local_range();
  for (unsigned int i = local_begin; i < local_end; ++i)
    to_fill(i) += value;
  to_fill.compress(dealii::VectorOperation::add);
}

// Verifies that the Updater takes ownership of the stamper.
TEST_F(IterationSourceUpdaterGaussSeidelTest, Constructor) {
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);

  EXPECT_EQ(mock_stamper_ptr_.use_count(), 2);
}

/*
 * ======== UpdateScatteringSource Tests =======================================
 */

// Verifies UpdateScatteringSource throws if RHS returns a null vector.
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateScatteringSourceBadRHS) {
  EXPECT_CALL(*mock_rhs_ptr_, GetVariableTermPtr(An<system::Index>(),_))
      .WillOnce(Return(nullptr));
  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);

  EXPECT_ANY_THROW(test_updater.UpdateScatteringSource(test_system_, 0, 0));
}

// Verify trying to update a group that has no moment returns an error
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateScatteringSourceBadMoment) {
  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);

  EXPECT_ANY_THROW(test_updater.UpdateScatteringSource(test_system_, 10, 0));
}

// Verifies operation of the UpdateScatteringSource function
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateScatteringSourceTestMPI) {
  system::GroupNumber group = test_helpers::RandomDouble(1, 3);
  system::AngleIndex angle = test_helpers::RandomDouble(0, 10);
  system::Index index = {group, angle};
  // Fill source vector with the value 2
  StampMPIVector(*source_vector_ptr_, 3);
  StampMPIVector(expected_vector_, group);


  /* Call expectations, expect to retrieve the scattering term vector from RHS
   * and then stamp it. We invoke the StampMPIVector function, which STAMPS a
   * vector. We make sure that the original value of 3, filled above, was zerod
   * out and replaced by the random group number.
   */
  EXPECT_CALL(*mock_rhs_ptr_, GetVariableTermPtr(index,
                                                 VariableTerms::kScatteringSource))
      .WillOnce(DoDefault());

  EXPECT_CALL(*moments_obs_ptr_, moments())
      .WillOnce(ReturnRef(current_iteration_moments_));

  EXPECT_CALL(*mock_stamper_ptr_,
      StampScatteringSource(Ref(*source_vector_ptr_),
                            group,
                            Ref(current_iteration_moments_)))
      .WillOnce(WithArgs<0,1>(Invoke(StampMPIVector)));

  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);
  // Tested call
  test_updater.UpdateScatteringSource(test_system_, group, angle);

  EXPECT_TRUE(bart::test_helpers::CompareMPIVectors(*source_vector_ptr_,
                                               expected_vector_));
}

/*
 * ======== UpdateFissionSource Tests ==========================================
 */

// Verifies UpdateFissionSource throws if RHS returns a null vector.
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateFissionSourceBadRHS) {
  EXPECT_CALL(*mock_rhs_ptr_, GetVariableTermPtr(An<system::Index>(),_))
      .WillOnce(Return(nullptr));
  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  test_system_.k_effective = 1.0;
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);

  EXPECT_ANY_THROW(test_updater.UpdateFissionSource(test_system_, 0, 0));
}

// Verify trying to update a group that has no moment returns an error
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateFissionSourceBadMoment) {
  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);

  EXPECT_ANY_THROW(test_updater.UpdateFissionSource(test_system_, 10, 0));
}

// Verify a bad keffective value (0, negative, or nullopt) returns an error
TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateFissionSourceBadKeff) {
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);

  // k_effective is still nullopt
  EXPECT_ANY_THROW(test_updater.UpdateFissionSource(test_system_, 0, 0));

  test_system_.k_effective = 0;
  EXPECT_ANY_THROW(test_updater.UpdateFissionSource(test_system_, 0, 0));

  test_system_.k_effective = -1;
  EXPECT_ANY_THROW(test_updater.UpdateFissionSource(test_system_, 0, 0));
}

TEST_F(IterationSourceUpdaterGaussSeidelTest, UpdateFissionSourceTestMPI) {
  system::GroupNumber group = test_helpers::RandomDouble(1, 3);
  system::AngleIndex angle = test_helpers::RandomDouble(0, 10);
  system::Index index = {group, angle};
  // Fill source vector with the value 3
  StampMPIVector(*source_vector_ptr_, 3);
  // Expected value identical to the group, which is [1, 3)
  StampMPIVector(expected_vector_, group);
  double k_effective = 1.05;
  test_system_.k_effective = k_effective;


  /* Call expectations, expect to retrieve the fission term vector from RHS
   * and then stamp it. We invoke the StampMPIVector function, which STAMPS a
   * vector. We make sure that the original value of 3, filled above, was zerod
   * out and replaced by the random group number.
   */
  EXPECT_CALL(*mock_rhs_ptr_, GetVariableTermPtr(index,
                                                 VariableTerms::kFissionSource))
      .WillOnce(DoDefault());

  std::array<int, 3> moment_index{group, 0, 0};
  EXPECT_CALL(*moments_obs_ptr_, GetMoment(moment_index))
      .WillOnce(ReturnRef(current_iteration_moments_[{group, 0, 0}]));
  EXPECT_CALL(*moments_obs_ptr_, moments())
      .WillOnce(ReturnRef(current_iteration_moments_));

  EXPECT_CALL(*mock_stamper_ptr_,
              StampFissionSource(Ref(*source_vector_ptr_),
                                 group,
                                 k_effective,
                                 Ref(current_iteration_moments_[{group, 0, 0}]),
                                 Ref(current_iteration_moments_)))
      .WillOnce(WithArgs<0,1>(Invoke(StampMPIVector)));

  // Final Set up
  test_system_.right_hand_side_ptr_ = std::move(mock_rhs_ptr_);
  CFEMSourceUpdater test_updater(mock_stamper_ptr_);
  // Tested call
  test_updater.UpdateFissionSource(test_system_, group, angle);

  EXPECT_TRUE(bart::test_helpers::CompareMPIVectors(*source_vector_ptr_,
                                               expected_vector_));
}

} // namespace