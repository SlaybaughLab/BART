#include "iteration/outer/outer_fixed_source_iteration.hpp"

#include "convergence/tests/iteration_completion_checker_mock.hpp"
#include "instrumentation/tests/instrument_mock.h"
#include "iteration/group/tests/group_solve_iteration_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"


namespace  {

using namespace bart;
using ::testing::A, ::testing::Ref, ::testing::_;

/* This fixture tests the operation of the OuterFixedSourceIteration class. This is a mediator class so the tests will verify
 * proper mediation between the dependencies and exposure of data to instruments.
 *
 * At the completion of SetUp() the test_iterator pointee object is set up with mock dependencies accessible via
 * the provided observation pointers. Mock instrumentation is accessible via shared pointers.
 * */
class IterationOuterFixedSourceIterationTest : public ::testing::Test {
 public:
  using TestIterationType = iteration::outer::OuterFixedSourceIteration;
  // Dependency types
  using GroupIteratorType = iteration::group::GroupSolveIterationMock;
  using ConvergenceCheckerType = typename convergence::IterationCompletionCheckerMock<double>;
  using ConvergenceInstrumentType = instrumentation::InstrumentMock<convergence::Status>;
  using StatusInstrumentType = instrumentation::InstrumentMock<std::string>;

  // Test object
  std::unique_ptr<TestIterationType> test_iterator{ nullptr };

  // Mock instruments
  std::shared_ptr<ConvergenceInstrumentType> convergence_status_instrument_ptr_ { std::make_shared<ConvergenceInstrumentType>() };
  std::shared_ptr<StatusInstrumentType> status_instrument_ptr_{ std::make_shared<StatusInstrumentType>() };

  // Observation pointers
  GroupIteratorType* group_iterator_mock_obs_ptr_;
  ConvergenceCheckerType* convergence_checker_mock_obs_ptr_;

  // Test parameters and objects
  bart::system::System test_system;

  void SetUp() override;
};

void IterationOuterFixedSourceIterationTest::SetUp() {
  auto group_iterator_ptr = std::make_unique<GroupIteratorType>();
  group_iterator_mock_obs_ptr_ = group_iterator_ptr.get();
  auto convergence_checker_ptr = std::make_unique<ConvergenceCheckerType>();
  convergence_checker_mock_obs_ptr_ = convergence_checker_ptr.get();

  test_iterator = std::make_unique<TestIterationType>(std::move(group_iterator_ptr),
                                                      std::move(convergence_checker_ptr));

  // Data ports
  using ConvergenceStatusPort = iteration::outer::data_names::ConvergenceStatusPort;
  using StatusPort = iteration::outer::data_names::StatusPort;

  instrumentation::GetPort<ConvergenceStatusPort>(*test_iterator).AddInstrument(convergence_status_instrument_ptr_);
  instrumentation::GetPort<StatusPort>(*test_iterator).AddInstrument(status_instrument_ptr_);
}

/* Constructor should not throw when valid dependency pointers are passed */
TEST_F(IterationOuterFixedSourceIterationTest, Constructor) {
  EXPECT_NO_THROW({
    TestIterationType test_iterator(std::make_unique<GroupIteratorType>(),std::make_unique<ConvergenceCheckerType>());
  });
}

/* Call to Iterate() should mediate the dependencies properly */
TEST_F(IterationOuterFixedSourceIterationTest, Iterate) {
  EXPECT_CALL(*convergence_checker_mock_obs_ptr_, Reset());
  EXPECT_CALL(*group_iterator_mock_obs_ptr_, Iterate(Ref(test_system))).Times(1);
  EXPECT_CALL(*convergence_status_instrument_ptr_, Read(_));
  EXPECT_CALL(*status_instrument_ptr_, Read(_));
  test_iterator->IterateToConvergence(this->test_system);
}

/* Call to UpdateSystem() should do nothing */
TEST_F(IterationOuterFixedSourceIterationTest, Update) {
  test_iterator->UpdateSystem(this->test_system, test_helpers::RandomInt(0, 10), test_helpers::RandomInt(0, 10));
}

} // namespace
