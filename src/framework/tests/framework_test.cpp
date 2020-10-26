#include "framework/framework.hpp"

#include <iostream>

#include <deal.II/base/mpi.h>

#include "iteration/outer/tests/outer_iteration_mock.hpp"
#include "iteration/initializer/tests/initializer_mock.h"
#include "results/tests/output_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::Ref, ::testing::_, ::testing::Return;

class FrameworkTest : public ::testing::Test {
 public:
  using Framework = framework::Framework;
  using Initializer = iteration::initializer::InitializerMock;
  using OuterIterator = iteration::outer::OuterIterationMock;
  using ResultsOutput = results::OutputMock;

  std::unique_ptr<Framework> test_framework_;

  // Observation Pointers
  system::System* system_obs_ptr_;
  Initializer* initializer_obs_ptr_;
  OuterIterator* outer_iterator_obs_ptr_;
  ResultsOutput* results_output_obs_ptr_;

  void SetUp() override;
};

void FrameworkTest::SetUp() {
  auto system_ptr = std::make_unique<system::System>();
  system_obs_ptr_ = system_ptr.get();
  auto initializer_ptr = std::make_unique<Initializer>();
  initializer_obs_ptr_ = initializer_ptr.get();
  auto outer_iterator_ptr = std::make_unique<OuterIterator>();
  outer_iterator_obs_ptr_ = outer_iterator_ptr.get();
  auto results_output_ptr = std::make_unique<ResultsOutput>();
  results_output_obs_ptr_ = results_output_ptr.get();

  test_framework_ = std::make_unique<Framework>(
      std::move(system_ptr),
      std::move(initializer_ptr),
      std::move(outer_iterator_ptr),
      std::move(results_output_ptr)
      );
}

TEST_F(FrameworkTest, Constructor) {
  EXPECT_NE(test_framework_->system(), nullptr);
  EXPECT_NE(test_framework_->initializer_ptr(), nullptr);
  EXPECT_NE(test_framework_->outer_iterator_ptr(), nullptr);
  EXPECT_NE(test_framework_->results_output_ptr(), nullptr);
}

TEST_F(FrameworkTest, ConstructorThrows) {
  for (int i = 0; i < 3; ++i) {
    auto system_ptr = (i == 0) ? nullptr :
        std::make_unique<system::System>();
    auto initializer_ptr = (i == 1) ? nullptr :
        std::make_unique<Initializer>();
    auto outer_iterator_ptr = (i == 2) ? nullptr :
        std::make_unique<OuterIterator>();
    auto results_output_ptr = std::make_unique<ResultsOutput>();
    EXPECT_ANY_THROW({
                       Framework test_framework(
                           std::move(system_ptr),
                           std::move(initializer_ptr),
                           std::move(outer_iterator_ptr),
                           std::move(results_output_ptr));
    });
  }
}

TEST_F(FrameworkTest, SolveSystem) {
  auto& system = *system_obs_ptr_;

  EXPECT_CALL(*initializer_obs_ptr_, Initialize(Ref(system)));
  EXPECT_CALL(*outer_iterator_obs_ptr_, IterateToConvergence(Ref(system)));

  test_framework_->SolveSystem();
}

TEST_F(FrameworkTest, OutputResultsNoOutputter) {
  auto system_ptr = std::make_unique<system::System>();
  auto initializer_ptr = std::make_unique<Initializer>();
  auto outer_iterator_ptr = std::make_unique<OuterIterator>();

  Framework test_framework(std::move(system_ptr),
                           std::move(initializer_ptr),
                           std::move(outer_iterator_ptr));

  std::ostringstream output_stream;

  EXPECT_ANY_THROW(test_framework.OutputResults(output_stream));
}

TEST_F(FrameworkTest, OutputResultsBadStream) {
  std::ostringstream output_stream;
  output_stream.setstate(std::ios_base::badbit);

  EXPECT_ANY_THROW(test_framework_->OutputResults(output_stream));
}

TEST_F(FrameworkTest, OutputResults) {
  std::ostringstream output_stream;
  auto& system = *system_obs_ptr_;
  auto& results_output_mock = *results_output_obs_ptr_;

  EXPECT_CALL(results_output_mock, AddData(Ref(system)));
  EXPECT_CALL(*results_output_obs_ptr_, WriteData(Ref(output_stream)));

  test_framework_->OutputResults(output_stream);
}

TEST_F(FrameworkTest, OutputMasterFileMPI) {
  std::ostringstream output_stream;
  const int process_id = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::vector<std::string> filenames{"mock_1", "mock_2"};

  if (process_id == 0) {
    EXPECT_CALL(*results_output_obs_ptr_, WriteMasterFile(Ref(output_stream),
                                                          filenames));
  } else {
    EXPECT_CALL(*results_output_obs_ptr_, WriteMasterFile(_,_))
        .Times(0);
  }

  test_framework_->OutputMasterFile(output_stream, filenames, process_id);
}

} // namespace