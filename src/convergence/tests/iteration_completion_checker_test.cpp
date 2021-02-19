#include "convergence/iteration_completion_checker.hpp"
#include "convergence/tests/convergence_checker_mock.hpp"
#include "convergence/status.hpp"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::Expectation, ::testing::NiceMock, ::testing::Return;

template <typename CompareT>
class ConvergenceIterationCompletionCheckerTest : public ::testing::Test {
 public:
  using ConvergenceCheckerMock = NiceMock<convergence::ConvergenceCheckerMock<CompareT>>;
  using IterationCompletionChecker = convergence::IterationCompletionChecker<CompareT>;

  std::unique_ptr<IterationCompletionChecker> test_completion_checker_ptr_;

  ConvergenceCheckerMock* convergence_checker_mock_obs_ptr_;
  auto SetUp() -> void override;
};

template <typename CompareT>
auto ConvergenceIterationCompletionCheckerTest<CompareT>::SetUp() -> void {
  auto convergence_checker_mock_ptr = std::make_unique<ConvergenceCheckerMock>();
  convergence_checker_mock_obs_ptr_ = convergence_checker_mock_ptr.get();
  test_completion_checker_ptr_ = std::make_unique<IterationCompletionChecker>(std::move(convergence_checker_mock_ptr));
}

using TestTypes = ::testing::Types<double>;
TYPED_TEST_SUITE(ConvergenceIterationCompletionCheckerTest, TestTypes);

// Getters should return the correct values, all convergence status values should be set to whatever
// the default values for a status are.
TYPED_TEST(ConvergenceIterationCompletionCheckerTest, Getters) {
  auto& test_completion_checker = *this->test_completion_checker_ptr_;
  const convergence::Status default_status;
  EXPECT_EQ(test_completion_checker.convergence_checker_ptr(), this->convergence_checker_mock_obs_ptr_);
  EXPECT_EQ(test_completion_checker.convergence_status(), default_status);
  EXPECT_EQ(test_completion_checker.max_iterations(), default_status.max_iterations);
  EXPECT_EQ(test_completion_checker.iteration(), default_status.iteration_number);
  EXPECT_EQ(test_completion_checker.convergence_is_complete(), default_status.is_complete);
}

// Setters should set the correct values
TYPED_TEST(ConvergenceIterationCompletionCheckerTest, Setters) {
  auto& test_completion_checker = *this->test_completion_checker_ptr_;
  const convergence::Status default_status;
  EXPECT_EQ(test_completion_checker.max_iterations(), default_status.max_iterations);
  const int new_max_iterations{ test_helpers::RandomInt(100, 200) };
  test_completion_checker.SetMaxIterations(new_max_iterations);
  EXPECT_EQ(test_completion_checker.max_iterations(), new_max_iterations);
  const int new_iteration{ test_helpers::RandomInt(0, 50) };
  test_completion_checker.SetIteration(new_iteration);
  EXPECT_EQ(test_completion_checker.iteration(), new_iteration);
}

// Setter max iterations to 0 or negative should throw
TYPED_TEST(ConvergenceIterationCompletionCheckerTest, SetMaxIterationsBadValues) {
  for (const int bad_max_iteration_value : {-10, -20, 0}) {
    EXPECT_ANY_THROW(this->test_completion_checker_ptr_->SetMaxIterations(bad_max_iteration_value));
  }
}

// Setter iterations to negative should throw
TYPED_TEST(ConvergenceIterationCompletionCheckerTest, SetIterationsBadValues) {
  for (const int bad_iteration_value : {-10, -20}) {
    EXPECT_ANY_THROW(this->test_completion_checker_ptr_->SetIteration(bad_iteration_value));
  }
}

TYPED_TEST(ConvergenceIterationCompletionCheckerTest, BadConvergenceAfterNSteps) {
  const int n_good_steps{test_helpers::RandomInt(1, 20) };
  auto good_previous_values {test_helpers::RandomVector(n_good_steps, -100, 100) };
  auto good_current_values {test_helpers::RandomVector(n_good_steps, -100, 100) };
  double bad_previous_value{test_helpers::RandomDouble(100, 200) };
  double bad_current_value{bad_previous_value + 1 };
  const double bad_delta{test_helpers::RandomDouble(-100, 100) };
  const int bad_index{ test_helpers::RandomInt(0, 100) };

  for (auto n = 0; n < n_good_steps; ++n) {
    Expectation good_step = EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_,
                                        IsConverged(good_current_values.at(n), good_previous_values.at(n)))
        .WillOnce(Return(true));
    EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, delta())
        .After(good_step)
        .WillOnce(Return(test_helpers::RandomDouble(-100, 100)));
  }

  Expectation bad_convergence = EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_,
                                            IsConverged(bad_current_value, bad_previous_value))
      .WillOnce(Return(false));
  EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, delta())
      .After(bad_convergence)
      .WillOnce(Return(bad_delta));
  EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, failed_index())
      .After(bad_convergence)
      .WillOnce(Return(bad_index));

  convergence::Status result;
  const convergence::Status expected_status{ .iteration_number = n_good_steps + 1,
                                             .max_iterations = result.max_iterations, // Ensures set to the default
                                             .is_complete = false,
                                             .failed_index = bad_index,
                                             .delta = bad_delta };
  for (auto n = 0; n < n_good_steps; ++n) {
    result = this->test_completion_checker_ptr_->ConvergenceStatus(good_current_values.at(n),
                                                                   good_previous_values.at(n));
  }
  result = this->test_completion_checker_ptr_->ConvergenceStatus(bad_current_value, bad_previous_value);
  EXPECT_EQ(result, expected_status);
}

TYPED_TEST(ConvergenceIterationCompletionCheckerTest, GoodConvergenceAfterNSteps) {
  const int n_bad_steps{test_helpers::RandomInt(1, 20) };
  auto bad_previous_values {test_helpers::RandomVector(n_bad_steps, -100, 100) };
  auto bad_current_values {test_helpers::RandomVector(n_bad_steps, -100, 100) };
  double good_previous_value{test_helpers::RandomDouble(100, 200) };
  double good_current_value{good_previous_value + 1 };
  const double good_delta{test_helpers::RandomDouble(-100, 100) };

  for (auto n = 0; n < n_bad_steps; ++n) {
    Expectation bad_convergence = EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_,
                                              IsConverged(bad_current_values.at(n), bad_previous_values.at(n)))
        .WillOnce(Return(false));
    EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, delta())
        .After(bad_convergence)
        .WillOnce(Return(test_helpers::RandomDouble(-100, 100)));
    EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, failed_index())
        .After(bad_convergence)
        .WillOnce(Return(test_helpers::RandomInt(0, 100)));
  }

  Expectation good_convergence = EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_,
                                            IsConverged(good_current_value, good_previous_value))
      .WillOnce(Return(true));
  EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, delta())
      .After(good_convergence)
      .WillOnce(Return(good_delta));

  convergence::Status result;
  const convergence::Status expected_status{ .iteration_number = n_bad_steps + 1,
      .max_iterations = result.max_iterations, // Ensures set to the default
      .is_complete = true,
      .failed_index = std::nullopt,
      .delta = good_delta };
  for (auto n = 0; n < n_bad_steps; ++n) {
    result = this->test_completion_checker_ptr_->ConvergenceStatus(bad_current_values.at(n),
                                                                   bad_previous_values.at(n));
  }
  result = this->test_completion_checker_ptr_->ConvergenceStatus(good_current_value, good_previous_value);
  EXPECT_EQ(result, expected_status);
}

TYPED_TEST(ConvergenceIterationCompletionCheckerTest, ReportConvergenceAfterMaxSteps) {
  const int max_iterations{ test_helpers::RandomInt(1, 10) };
  auto bad_previous_values { test_helpers::RandomVector(max_iterations, -100, 100) };
  auto bad_current_values { test_helpers::RandomVector(max_iterations, -100, 100) };
  std::vector<double> delta;
  std::vector<int> failed_indices;

  for (auto n = 0; n < max_iterations; ++n) {
    Expectation step = EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_,
                                   IsConverged(bad_current_values.at(n), bad_previous_values.at(n)))
        .WillOnce(Return(false));
    const double delta_value{ test_helpers::RandomDouble(-100, 100) };
    const int failed_index{ test_helpers::RandomInt(0, 100) };
    delta.push_back(delta_value);
    failed_indices.push_back(failed_index);
    EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, delta())
        .After(step).WillOnce(Return(delta_value));
    EXPECT_CALL(*this->convergence_checker_mock_obs_ptr_, failed_index())
        .After(step).WillOnce(Return(failed_index));
  }
  this->test_completion_checker_ptr_->SetMaxIterations(max_iterations);
  convergence::Status result;

  for (auto n = 0; n < max_iterations; ++n) {
    result = this->test_completion_checker_ptr_->ConvergenceStatus(bad_current_values.at(n),
                                                                   bad_previous_values.at(n));
  }

  EXPECT_EQ(result.max_iterations, max_iterations);
  EXPECT_EQ(result.iteration_number, max_iterations);
  EXPECT_TRUE(result.is_complete);
  EXPECT_NE(std::find(delta.cbegin(), delta.cend(), result.delta), delta.end());
  EXPECT_TRUE(result.failed_index.has_value());
  EXPECT_NE(std::find(failed_indices.cbegin(), failed_indices.cend(), result.failed_index), failed_indices.end());
}





} // namespace
