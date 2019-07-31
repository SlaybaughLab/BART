#include "iteration/group/group_source_iteration.h"

#include <functional>
#include <memory>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_full_matrix.h>

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
using ::testing::AtLeast;
using ::testing::ExpectationSet;
using ::testing::Return, ::testing::Pointee, ::testing::Ref;
using ::testing::ReturnRef;
using ::testing::Sequence, ::testing::_;
using ::testing::InvokeWithoutArgs;
using ::testing::Unused;

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

  virtual ~IterationGroupSourceIterationTest() = default;

  // Test object
  std::unique_ptr<TestGroupIterator> test_iterator_ptr_;

  // Mock dependency objects
  std::unique_ptr<GroupSolver> single_group_solver_ptr_;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr_;
  std::shared_ptr<GroupSolution> group_solution_ptr_;
  std::shared_ptr<SourceUpdater> source_updater_ptr_;

  // Supporting objects
  system::System test_system;

  // Observing pointers
  GroupSolver* single_group_obs_ptr_ = nullptr;
  ConvergenceChecker* convergence_checker_obs_ptr_ = nullptr;
  MomentCalculator* moment_calculator_obs_ptr_ = nullptr;
  SourceUpdater* source_updater_obs_ptr_ = nullptr;
  Moments* moments_obs_ptr_ = nullptr;

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
  source_updater_ptr_ = std::make_shared<SourceUpdater>();
  source_updater_obs_ptr_ = source_updater_ptr_.get();

  test_system.current_moments = std::make_unique<Moments>();
  moments_obs_ptr_ = dynamic_cast<Moments*>(test_system.current_moments.get());

  test_iterator_ptr_ = std::make_unique<TestGroupIterator>(
      std::move(single_group_solver_ptr_),
      std::move(convergence_checker_ptr_),
      std::move(moment_calculator_ptr_),
      group_solution_ptr_,
      source_updater_ptr_
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

template <typename DimensionWrapper>
class IterationGroupSourceSystemSolvingTest :
    public IterationGroupSourceIterationTest<DimensionWrapper> {
 public:
  virtual ~IterationGroupSourceSystemSolvingTest() = default;
  IterationGroupSourceSystemSolvingTest()
      : L_(4,4), U_(4,4), b_(MPI_COMM_WORLD, 4, 4),
        true_scalar_flux_(4),
        solver_(solver_control_, MPI_COMM_WORLD)
        {}
  // Test parameters
  static constexpr int total_groups = 2;
  static constexpr int total_angles = 2;
  static constexpr int max_harmonic_l = 1;

  // Objects for solver
  dealii::PETScWrappers::FullMatrix L_;
  dealii::PETScWrappers::FullMatrix U_;
  dealii::PETScWrappers::MPI::Vector b_;
  system::moments::MomentVector true_scalar_flux_;

  // Solver
  dealii::SolverControl solver_control_;
  dealii::PETScWrappers::SolverGMRES solver_;

  // Group solutions
  std::array<dealii::PETScWrappers::MPI::Vector, total_groups> group_solutions_;
  std::array<dealii::PETScWrappers::MPI::Vector, total_groups> group_rhs_;

  int iterations = 0;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationGroupSourceSystemSolvingTest<DimensionWrapper>::SetUp() {
  IterationGroupSourceIterationTest<DimensionWrapper>::SetUp();
  // SYSTEM SETUP ==============================================================
  // Matrix and vectors for system to be solved
  const std::array<double, 16> l_entries{  2.69110073,  0.        ,  0.        ,  0.        ,
                                           -2.27591718, 13.46409688,  0.        ,  0.        ,
                                           1.50909047, -3.12175911,  4.01996062,  0.        ,
                                           3.15786087, -4.17175472, -1.10106215,  9.55326528};
  const std::array<double, 16> u_entries{ 0.        , -2.27591718,  1.50909047,  3.15786087,
                                          0.        ,  0.        , -3.12175911, -4.17175472,
                                          0.        ,  0.        ,  0.        , -1.10106215,
                                          0.        ,  0.        ,  0.        ,  0.        };
  const std::array<double, 4> b_entries{ 4.53187178, -4.15513557, 6.79740634,
                                         2.43451537};
  // Solution
  const std::array<double, 4> true_scalar_flux_entries{-0.03980115, 0.42018347, 2.2260879,
                                                       0.70804747};

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      L_.set(i,j, l_entries.at(4*i + j));
      U_.set(i,j, u_entries.at(4*i + j));
    }
    b_[i] = b_entries.at(i);
    true_scalar_flux_[i] = true_scalar_flux_entries.at(i);
  }
  L_.compress(dealii::VectorOperation::insert);
  U_.compress(dealii::VectorOperation::insert);
  b_.compress(dealii::VectorOperation::insert);

  for (int group = 0; group < this->total_groups; ++group) {
    dealii::PETScWrappers::MPI::Vector group_solution(MPI_COMM_WORLD, 4, 4),
        group_rhs_vector(MPI_COMM_WORLD, 4, 4);
    for (int i = 0; i < 4; ++i) {
      group_solution[i] = 1;
    }
    group_solution.compress(dealii::VectorOperation::insert);
    // Generate initial RHS, b - Ux_0
    U_.vmult_add(group_rhs_vector, group_solution);
    group_rhs_vector.add(-1, b_);
    group_rhs_vector *= -1;
    group_rhs_vector.compress(dealii::VectorOperation::add);

    group_solutions_[group] = std::move(group_solution);
    group_rhs_[group] = std::move(group_rhs_vector);

  }
}

ACTION_P(Solve, test_class) {
  dealii::PETScWrappers::PreconditionNone no_conditioner(test_class->L_);
  test_class->solver_.solve(test_class->L_,
                            test_class->group_solutions_.at(arg0),
                            test_class->group_rhs_.at(arg0),
                            no_conditioner);
}

ACTION_P(Update, test_class) {
  auto& rhs = test_class->group_rhs_.at(arg1);
  auto& solution = test_class->group_solutions_.at(arg1);
  rhs = 0;
  test_class->U_.vmult_add(rhs, solution);
  rhs.add(-1, test_class->b_);
  rhs *= -1;
}

ACTION_P(CalculatedScalarFlux, test_class) {
  system::moments::MomentVector return_vector(test_class->group_solutions_.at(arg1));
  return_vector *= (arg2 + arg3 + 1);
  return return_vector;
}

ACTION_P(ReturnConvergence, test_class) {
  convergence::Status status;
  ++test_class->iterations;
  auto diff = arg0;
  diff.add(-1, arg1);
  double err = diff.l1_norm();
  if (err < 1e-6 || test_class->iterations > 1000)
    status.is_complete = true;
  return status;
}

ACTION_P(ResetIterations, test_class) {
  test_class->iterations = 0;
}

TYPED_TEST_CASE(IterationGroupSourceSystemSolvingTest, bart::testing::AllDimensions);

TYPED_TEST(IterationGroupSourceSystemSolvingTest, Iterate) {
  // This is the mock map to hold system current_moments
  system::moments::MomentsMap current_moments;
  for (int group = 0; group < this->total_groups; ++group) {
    for (int l = 0; l <= this->max_harmonic_l; ++l) {
      for (int m = -l; m <= l; ++m) {
        system::moments::MomentIndex index{group, l, m};
        current_moments.emplace(index, 4);
        EXPECT_CALL(*this->moments_obs_ptr_, BracketOp(index))
            .Times(AtLeast(1))
            .WillRepeatedly(ReturnRef(current_moments.at(index)));

        EXPECT_CALL(*this->moment_calculator_obs_ptr_, CalculateMoment(
            this->group_solution_ptr_.get(), group, l, m))
            .Times(AtLeast(1))
            .WillRepeatedly(CalculatedScalarFlux(this));
      }
    }
    EXPECT_CALL(*this->single_group_obs_ptr_, SolveGroup(
        group, Ref(this->test_system), Ref(*this->group_solution_ptr_)))
        .Times(AtLeast(1))
        .WillRepeatedly(Solve(this));

    for (int angle = 0; angle < this->total_angles; ++angle) {
      EXPECT_CALL(*this->source_updater_obs_ptr_, UpdateScatteringSource(
          Ref(this->test_system), group, angle))
          .Times(AtLeast(1))
          .WillRepeatedly(Update(this));
    }
  }

  EXPECT_CALL(*this->convergence_checker_obs_ptr_, Reset())
  .Times(AtLeast(1))
  .WillRepeatedly(ResetIterations(this));

  EXPECT_CALL(*this->convergence_checker_obs_ptr_, CheckFinalConvergence(_, _))
  .Times(AtLeast(1))
  .WillRepeatedly(ReturnConvergence(this));

  EXPECT_CALL(*this->moments_obs_ptr_, total_groups())
      .WillOnce(Return(this->total_groups));
  EXPECT_CALL(*this->group_solution_ptr_, total_angles())
      .WillOnce(Return(this->total_angles));
  EXPECT_CALL(*this->moments_obs_ptr_, max_harmonic_l())
      .WillRepeatedly(Return(this->max_harmonic_l));

  this->test_iterator_ptr_->Iterate(this->test_system);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(current_moments.at({0, 0, 0})[i],
                     this->true_scalar_flux_[i], 1e-6);
  }
}
} // namespace