#include "iteration/outer/outer_fixed_source_iteration.h"

#include "test_helpers/gmock_wrapper.h"

#include "iteration/group/tests/group_solve_iteration_mock.h"
#include "convergence/tests/final_checker_mock.h"
#include "convergence/reporter/tests/mpi_mock.h"

namespace  {

using namespace bart;

class IterationOuterDummyIterationTest : public ::testing::Test {
 public:
  using TestIterationType = iteration::outer::OuterFixedSourceIteration;
  // Dependency types
  using GroupIteratorType = iteration::group::GroupSolveIterationMock;
  using ConvergenceCheckerType = typename convergence::FinalCheckerMock<double>;
  using ReporterType = convergence::reporter::MpiMock;

  // Test object
  std::unique_ptr<TestIterationType> test_iterator = nullptr;

  // Dependencies
  std::shared_ptr<ReporterType> reporter_mock_ptr_;

  // Observation pointers
  GroupIteratorType* group_iterator_mock_obs_ptr_;
  ConvergenceCheckerType* convergence_checker_mock_obs_ptr_;

  void SetUp() override;
};

void IterationOuterDummyIterationTest::SetUp() {
  reporter_mock_ptr_ = std::make_shared<ReporterType>();
  auto group_iterator_ptr = std::make_unique<GroupIteratorType>();
  group_iterator_mock_obs_ptr_ = group_iterator_ptr.get();
  auto convergence_checker_ptr = std::make_unique<ConvergenceCheckerType>();
  convergence_checker_mock_obs_ptr_ = convergence_checker_ptr.get();

  test_iterator = std::make_unique<TestIterationType>(
      std::move(group_iterator_ptr),
      std::move(convergence_checker_ptr),
      reporter_mock_ptr_);
}

TEST_F(IterationOuterDummyIterationTest, Constructor) {
  EXPECT_NO_THROW({
    TestIterationType test_iterator(
        std::make_unique<GroupIteratorType>(),
        std::make_unique<ConvergenceCheckerType>(),
        reporter_mock_ptr_
        );
  });
}

} // namespace
