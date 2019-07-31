#include "iteration/outer/outer_power_iteration.h"

#include <memory>

#include "iteration/group/tests/group_solve_iteration_mock.h"
#include "eigenvalue/k_effective/tests/k_effective_updater_mock.h"
#include "convergence/reporter/tests/mpi_mock.h"
#include "convergence/tests/final_checker_mock.h"
#include "iteration/updater/tests/source_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "system/system.h"

namespace  {

using namespace bart;

using ::testing::A, ::testing::AtLeast, ::testing::Expectation;
using ::testing::Ref, ::testing::Return, ::testing::Sequence, ::testing::_;

class IterationOuterPowerIterationTest : public ::testing::Test {
 protected:
  using GroupIterator = iteration::group::GroupSolveIterationMock;
  using ConvergenceChecker = convergence::FinalCheckerMock<double>;
  using K_EffectiveUpdater = eigenvalue::k_effective::K_EffectiveUpdaterMock;
  using OuterPowerIteration = iteration::outer::OuterPowerIteration;
  using SourceUpdater = iteration::updater::SourceUpdaterMock;
  using Reporter = convergence::reporter::MpiMock;

  std::unique_ptr<OuterPowerIteration> test_iterator;

  // Dependencies
  std::shared_ptr<SourceUpdater> source_updater_ptr_;
  std::shared_ptr<Reporter> reporter_ptr_;

  // Supporting objects
  system::System test_system;

  // Observation pointers
  GroupIterator* group_iterator_obs_ptr_;
  ConvergenceChecker* convergence_checker_obs_ptr_;
  K_EffectiveUpdater* k_effective_updater_obs_ptr_;

  // Test parameters
  const int total_groups = 2;
  const int total_angles = 3;
  static constexpr int iterations_ = 4;

  void SetUp() override;
};

void IterationOuterPowerIterationTest::SetUp() {
  // Dependencies
  source_updater_ptr_ = std::make_shared<SourceUpdater>();
  auto group_iterator_ptr = std::make_unique<GroupIterator>();
  group_iterator_obs_ptr_ = group_iterator_ptr.get();
  auto convergenge_checker_ptr = std::make_unique<ConvergenceChecker>();
  convergence_checker_obs_ptr_ = convergenge_checker_ptr.get();
  auto k_effective_updater_ptr = std::make_unique<K_EffectiveUpdater>();
  k_effective_updater_obs_ptr_ = k_effective_updater_ptr.get();
  reporter_ptr_ = std::make_shared<Reporter>();

  // Set up system
  test_system.total_angles = total_angles;
  test_system.total_groups = total_groups;

  // Construct test object
  test_iterator = std::make_unique<OuterPowerIteration>(
      std::move(group_iterator_ptr),
      std::move(convergenge_checker_ptr),
      std::move(k_effective_updater_ptr),
      source_updater_ptr_,
      reporter_ptr_
      );
}

TEST_F(IterationOuterPowerIterationTest, Constructor) {
  EXPECT_NE(this->test_iterator, nullptr);
  EXPECT_NE(this->test_iterator->group_iterator_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->source_updater_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->convergence_checker_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->k_effective_updater_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->reporter_ptr(), nullptr);
  EXPECT_EQ(this->source_updater_ptr_.use_count(), 2);
}

TEST_F(IterationOuterPowerIterationTest, ConstructorErrors) {

  for (int i = 0; i < 4; ++i) {
    auto convergence_checker_ptr = (i == 0) ? nullptr :
        std::make_unique<convergence::FinalCheckerMock<double>>();
    auto k_effective_updater_ptr = (i == 1) ? nullptr :
        std::make_unique<eigenvalue::k_effective::K_EffectiveUpdaterMock>();
    auto source_updater_ptr = (i == 2) ? nullptr : this->source_updater_ptr_;
    auto group_iterator_ptr = (i == 3) ? nullptr :
        std::make_unique<iteration::group::GroupSolveIterationMock>();

    EXPECT_ANY_THROW({
                       iteration::outer::OuterPowerIteration test_iterator(
                           std::move(group_iterator_ptr),
                           std::move(convergence_checker_ptr),
                           std::move(k_effective_updater_ptr),
                           source_updater_ptr);
                     });
  }
}

TEST_F(IterationOuterPowerIterationTest, IterateToConvergenceTest) {

  for (int group = 0; group < this->total_groups; ++group) {
    for (int angle = 0; angle < this->total_angles; ++angle) {
      EXPECT_CALL(*this->source_updater_ptr_, UpdateFissionSource(
          Ref(this->test_system),group, angle))
          .Times(this->iterations_);
    }
  }

  // K_Effective updater return values
  std::array<double, this->iterations_ + 1> k_effective_by_iteration;
  Sequence k_effective_calls;

  for (int i = 0; i < this->iterations_; ++i) {
    k_effective_by_iteration.at(i + 1) = i * 1.5;

    EXPECT_CALL(*this->k_effective_updater_obs_ptr_,
                CalculateK_Effective(Ref(this->test_system)))
                .InSequence(k_effective_calls)
                .WillOnce(Return(k_effective_by_iteration.at(i + 1)));

    convergence::Status convergence_status;
    convergence_status.is_complete = (i == (this->iterations_ - 1));

    EXPECT_CALL(*this->convergence_checker_obs_ptr_,
                CheckFinalConvergence(k_effective_by_iteration.at(i + 1),
                                      k_effective_by_iteration.at(i)))
            .WillOnce(Return(convergence_status));
  }

  EXPECT_CALL(*this->group_iterator_obs_ptr_, Iterate(Ref(this->test_system)))
      .Times(this->iterations_);

  EXPECT_CALL(*this->reporter_ptr_, Report(A<const convergence::Status&>()))
      .Times(this->iterations_);
  EXPECT_CALL(*this->reporter_ptr_, Report(A<const std::string&>()))
      .Times(AtLeast(this->iterations_));

  this->test_iterator->IterateToConvergence(this->test_system);
}



} // namespace