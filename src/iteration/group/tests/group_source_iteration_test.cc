#include "iteration/group/group_source_iteration.h"

#include <memory>

#include "iteration/updater/tests/source_updater_mock.h"
#include "quadrature/calculators/tests/spherical_harmonic_moments_mock.h"
#include "convergence/tests/final_checker_mock.h"
#include "solver/group/tests/single_group_solver_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/system.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return, ::testing::Pointee, ::testing::Ref;
using ::testing::ReturnRef;
using ::testing::Sequence, ::testing::_;

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
  using Moments = system::moments::SphericalHarmonicMock;

  // Test object
  std::unique_ptr<TestGroupIterator> test_iterator_ptr_;

  // Mock dependency objects
  std::unique_ptr<GroupSolver> single_group_solver_ptr_;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr_;
  std::shared_ptr<GroupSolution> group_solution_ptr_;
  std::unique_ptr<SourceUpdater> source_updater_ptr_;

  // Supporting objects
  system::System test_system;

  // Observing pointers
  GroupSolver* single_group_obs_ptr_ = nullptr;
  ConvergenceChecker* convergence_checker_obs_ptr_ = nullptr;
  MomentCalculator* moment_calculator_obs_ptr_ = nullptr;
  SourceUpdater* source_updater_obs_ptr_ = nullptr;
  Moments* moments_obs_ptr_ = nullptr;

  // Test parameters
  const int total_groups = 2;
  const int total_angles = 3;
  const std::array<int, 2> iterations_by_group{2,3};
  const int max_harmonic_l = 2;

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

  test_system.current_moments = std::make_unique<Moments>();
  moments_obs_ptr_ = dynamic_cast<Moments*>(test_system.current_moments.get());

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

TYPED_TEST(IterationGroupSourceIterationTest, Iterate) {

  // Objects to support testing
  /* Moment Vectors. These will be returned by the MomentCalculator, based on
   * group and iteration. It is important to ensure that the correct values
   * are being checked for convergence. All entries in each vector will be set
   * to a unique value, 10*group + iteration. */
  std::map<int, std::vector<system::moments::MomentVector>> calculated_moments;
  system::moments::MomentVector zero_moment(5);

  for (int group = 0; group < this->total_groups; ++group) {
    calculated_moments[group] = {};
    for (int it = 0; it < this->iterations_by_group[group]; ++it) {
      system::moments::MomentVector new_moment(5);
      new_moment = (group * 10 + it);
      calculated_moments.at(group).push_back(new_moment);
    }
  }

  /* Final Moment Vectors. These will be returned by the MomentCalculator at the
   * end of all calculations, to update all moments (The previous moments are
   * only for scalar flux) */
  std::map<int, system::moments::MomentsMap> moments_map;
  std::map<int, system::moments::MomentsMap> returned_moments;

  for (int group = 0; group < this->total_groups; ++group) {
    moments_map[group] = {};
    returned_moments[group] = {};
    for (int l = 0; l <= this->max_harmonic_l; ++l) {
      for (int m = -l; m <= this->max_harmonic_l; ++m) {
        system::moments::MomentVector new_moment(5), new_empty_moment(5);
        new_moment = (l * 10) + m;
        moments_map[group][{l,m}] = new_moment;
        returned_moments[group][{l,m}] = new_empty_moment;
      }
    }
  }

  EXPECT_CALL(*this->moments_obs_ptr_, total_groups())
      .WillOnce(Return(this->total_groups));
  EXPECT_CALL(*this->group_solution_ptr_, total_angles())
      .WillOnce(Return(this->total_angles));

  Sequence s;

  for (int group = 0; group < this->total_groups; ++group) {
    EXPECT_CALL(*this->single_group_obs_ptr_, SolveGroup(
        group,
        Ref(this->test_system),
        Ref(*this->group_solution_ptr_)))
        .Times(this->iterations_by_group[group]);

    for (int it = 0; it < this->iterations_by_group[group]; ++it) {

      EXPECT_CALL(*this->moment_calculator_obs_ptr_, CalculateMoment(
          this->group_solution_ptr_.get(), group, 0, 0))
          //_, group, 0, 0))
          .InSequence(s)
          .WillOnce(Return(calculated_moments.at(group).at(it)));

      convergence::Status status;

      if ((it + 1) == this->iterations_by_group[group])
        status.is_complete = true;

      status.iteration_number = it + 1;

      if (it == 0) {
        EXPECT_CALL(*this->convergence_checker_obs_ptr_, CheckFinalConvergence(
            calculated_moments.at(group).at(it),
            zero_moment))
            .InSequence(s)
            .WillOnce(Return(status));
      } else {
        EXPECT_CALL(*this->convergence_checker_obs_ptr_, CheckFinalConvergence(
            calculated_moments.at(group).at(it),
            calculated_moments.at(group).at(it - 1)))
            .InSequence(s)
            .WillOnce(Return(status));
      }
    }
    // Updates should occur iterations - 1 times (doesn't run if converged)
    for (int angle = 0; angle < this->total_angles; ++angle) {
      EXPECT_CALL(*this->source_updater_obs_ptr_, UpdateScatteringSource(
          Ref(this->test_system), group, angle))
          .Times(this->iterations_by_group[group] - 1);
    }
  }

  // Following all groups, expect to update moments
  EXPECT_CALL(*this->moments_obs_ptr_, max_harmonic_l())
      .Times(this->total_groups)
      .WillRepeatedly(Return(this->max_harmonic_l));

  for (int group = 0; group < this->total_groups; ++group) {
    for (int l = 0; l <= this->max_harmonic_l; ++l) {
      for (int m = -l; m <= l; ++m) {
        EXPECT_CALL(*this->moment_calculator_obs_ptr_, CalculateMoment(
            this->group_solution_ptr_.get(), group, l, m))
            .InSequence(s)
            .WillOnce(Return(moments_map.at(group).at({l, m})));
        system::moments::MomentIndex index{group,l,m};

        EXPECT_CALL(*this->moments_obs_ptr_, BracketOp(index))
            .WillOnce(ReturnRef(returned_moments.at(group).at({l,m})));
      }
    }
  }

  // Make sure correct vectors have been stored
  for (auto& group_moments : returned_moments) {
    auto& [group, group_moment_map] = group_moments;
    for (auto& moment : group_moment_map) {
      auto& [index, returned_moment] = moment;
      EXPECT_EQ(returned_moment, moments_map.at(group).at(index));
    }
  }

  this->test_iterator_ptr_->Iterate(this->test_system);
}

} // namespace