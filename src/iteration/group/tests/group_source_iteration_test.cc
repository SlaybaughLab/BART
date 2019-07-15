#include "iteration/group/group_source_iteration.h"

#include <memory>

#include "iteration/updater/tests/source_updater_mock.h"
#include "quadrature/calculators/tests/spherical_harmonic_moments_mock.h"
#include "convergence/tests/final_checker_mock.h"
#include "solver/group/tests/single_group_solver_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class IterationGroupSourceIterationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using TestGroupIterator = iteration::group::GroupSourceIteration<dim>;
  using GroupSolver = solver::group::SingleGroupSolverMock;
  using ConvergenceChecker = convergence::FinalCheckerMock<system::moments::MomentVector>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsMock<dim>;
  using GroupSolution = system::solution::MPIGroupAngularSolutionMock;
  using SourceUpdater = iteration::updater::SourceUpdaterMock;

  // Test object
  std::unique_ptr<TestGroupIterator> test_iterator_ptr_;

  // Mock dependency objects
  std::unique_ptr<GroupSolver> single_group_solver_ptr_;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr_;
  std::shared_ptr<GroupSolution> group_solution_ptr_;
  std::unique_ptr<SourceUpdater> source_updater_ptr_;

  // Observing pointers
  GroupSolver* single_group_obs_ptr_ = nullptr;
  ConvergenceChecker* convergence_checker_obs_ptr_ = nullptr;
  MomentCalculator* moment_calculator_obs_ptr_ = nullptr;
  SourceUpdater* source_updater_obs_ptr_ = nullptr;

  void SetUp() override;
};

TYPED_TEST_CASE(IterationGroupSourceIterationTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void IterationGroupSourceIterationTest<DimensionWrapper>::SetUp() {
  single_group_solver_ptr_ = std::make_unique<GroupSolver>();
  single_group_obs_ptr_ = single_group_solver_ptr_.get();
  convergence_checker_ptr_ = std::make_unique<ConvergenceChecker>();
  convergence_checker_obs_ptr_ = convergence_checker_ptr_.get();
  moment_calculator_ptr_ = std::make_unique<MomentCalculator>();
  moment_calculator_obs_ptr_ = moment_calculator_ptr_.get();
  group_solution_ptr_ = std::make_shared<GroupSolution>();
  source_updater_ptr_ = std::make_unique<SourceUpdater>();
  source_updater_obs_ptr_ = source_updater_ptr_.get();

  test_iterator_ptr_ = std::make_unique<TestGroupIterator>(
      std::move(single_group_solver_ptr_),
      std::move(convergence_checker_ptr_),
      std::move(moment_calculator_ptr_),
      group_solution_ptr_,
      std::move(source_updater_ptr_)
      );
}

TYPED_TEST(IterationGroupSourceIterationTest, Constructor) {
  using GroupSolver = solver::group::SingleGroupSolverMock;
  using ConvergenceChecker = convergence::FinalCheckerMock<system::moments::MomentVector>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsMock<this->dim>;
  using SourceUpdater = iteration::updater::SourceUpdaterMock;

  auto single_group_test_ptr = dynamic_cast<GroupSolver*>(
      this->test_iterator_ptr_->group_solver_ptr());
  auto convergence_checker_test_ptr = dynamic_cast<ConvergenceChecker*>(
      this->test_iterator_ptr_->convergence_checker_ptr());
  auto moment_calculator_test_ptr = dynamic_cast<MomentCalculator*>(
      this->test_iterator_ptr_->moment_calculator_ptr());
  auto source_updater_test_ptr = dynamic_cast<SourceUpdater*>(
      this->test_iterator_ptr_->source_updater_ptr());

  EXPECT_NE(nullptr, single_group_test_ptr);
  EXPECT_NE(nullptr, convergence_checker_test_ptr);
  EXPECT_NE(nullptr, moment_calculator_test_ptr);
  EXPECT_EQ(2, this->group_solution_ptr_.use_count());
  EXPECT_EQ(this->group_solution_ptr_.get(),
            this->test_iterator_ptr_->group_solution_ptr().get());
  EXPECT_NE(nullptr, source_updater_test_ptr);
}

} // namespace