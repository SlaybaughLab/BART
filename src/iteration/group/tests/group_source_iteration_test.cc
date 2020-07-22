#include "iteration/group/group_source_iteration.h"

#include <functional>
#include <memory>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_full_matrix.h>

#include "formulation/updater/tests/boundary_conditions_updater_mock.h"
#include "formulation/updater/tests/scattering_source_updater_mock.h"
#include "quadrature/calculators/tests/spherical_harmonic_moments_mock.h"
#include "convergence/tests/final_checker_mock.h"
#include "convergence/reporter/tests/mpi_mock.h"
#include "solver/group/tests/single_group_solver_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/system.h"
#include "system/solution/solution_types.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"

namespace  {

using namespace bart;
using ::testing::AtLeast;
using ::testing::ExpectationSet;
using ::testing::Return, ::testing::Pointee, ::testing::Ref;
using ::testing::ReturnRef;
using ::testing::Sequence, ::testing::_;
using ::testing::InvokeWithoutArgs;
using ::testing::Unused;
using ::testing::A;

template <typename DimensionWrapper>
class IterationGroupSourceIterationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using EnergyGroupToAngularSolutionPtrMap = bart::system::solution::EnergyGroupToAngularSolutionPtrMap;
  using TestGroupIterator = iteration::group::GroupSourceIteration<dim>;
  using GroupSolver = solver::group::SingleGroupSolverMock;
  using ConvergenceChecker = convergence::FinalCheckerMock<system::moments::MomentVector>;
  using MomentMapConvergenceChecker = convergence::FinalCheckerMock<const system::moments::MomentsMap>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsMock;
  using GroupSolution = system::solution::MPIGroupAngularSolutionMock;
  using BoundaryConditionsUpdater = formulation::updater::BoundaryConditionsUpdaterMock;
  using SourceUpdater = formulation::updater::ScatteringSourceUpdaterMock;
  using Moments = system::moments::SphericalHarmonicMock;
  using Reporter = convergence::reporter::MpiMock;

  virtual ~IterationGroupSourceIterationTest() = default;

  // Test object
  std::unique_ptr<TestGroupIterator> test_iterator_ptr_;

  // Mock dependency objects
  std::shared_ptr<GroupSolution> group_solution_ptr_;
  std::shared_ptr<BoundaryConditionsUpdater> boundary_conditions_updater_ptr_;
  std::shared_ptr<SourceUpdater> source_updater_ptr_;
  std::shared_ptr<Reporter> reporter_ptr_;

  // Supporting objects
  system::System test_system;
  EnergyGroupToAngularSolutionPtrMap energy_group_angular_solution_ptr_map_;

  // Observing pointers
  GroupSolver* single_group_obs_ptr_ = nullptr;
  ConvergenceChecker* convergence_checker_obs_ptr_ = nullptr;
  MomentCalculator* moment_calculator_obs_ptr_ = nullptr;
  MomentMapConvergenceChecker* moment_map_convergence_checker_obs_ptr_ = nullptr;
  Moments* moments_obs_ptr_ = nullptr;
  Moments* previous_moments_obs_ptr_ = nullptr;

  void SetUp() override;
};

TYPED_TEST_CASE(IterationGroupSourceIterationTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void IterationGroupSourceIterationTest<DimensionWrapper>::SetUp() {
  auto single_group_solver_ptr_ = std::make_unique<GroupSolver>();
  single_group_obs_ptr_ = single_group_solver_ptr_.get();
  auto convergence_checker_ptr_ = std::make_unique<ConvergenceChecker>();
  convergence_checker_obs_ptr_ = convergence_checker_ptr_.get();
  auto moment_calculator_ptr_ = std::make_unique<MomentCalculator>();
  moment_calculator_obs_ptr_ = moment_calculator_ptr_.get();
  auto moment_map_convergence_checker_ptr_ =
      std::make_unique<MomentMapConvergenceChecker>();
  moment_map_convergence_checker_obs_ptr_ = moment_map_convergence_checker_ptr_.get();

  group_solution_ptr_ = std::make_shared<GroupSolution>();
  boundary_conditions_updater_ptr_ = std::make_shared<BoundaryConditionsUpdater>();
  source_updater_ptr_ = std::make_shared<SourceUpdater>();
  reporter_ptr_ = std::make_shared<Reporter>();

  test_system.current_moments = std::make_unique<Moments>();
  moments_obs_ptr_ = dynamic_cast<Moments*>(test_system.current_moments.get());
  test_system.previous_moments = std::make_unique<Moments>();
  previous_moments_obs_ptr_ = dynamic_cast<Moments*>(test_system.previous_moments.get());


  test_iterator_ptr_ = std::make_unique<TestGroupIterator>(
      std::move(single_group_solver_ptr_),
      std::move(convergence_checker_ptr_),
      std::move(moment_calculator_ptr_),
      group_solution_ptr_,
      source_updater_ptr_,
      boundary_conditions_updater_ptr_,
      reporter_ptr_,
      std::move(moment_map_convergence_checker_ptr_));
}

TYPED_TEST(IterationGroupSourceIterationTest, Constructor) {
  using GroupSolver = solver::group::SingleGroupSolverMock;
  using ConvergenceChecker = convergence::FinalCheckerMock<system::moments::MomentVector>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsMock;
  using MomentMapConvergenceChecker = convergence::FinalCheckerMock<const bart::system::moments::MomentsMap>;
  using SourceUpdater = formulation::updater::ScatteringSourceUpdaterMock;
  using BoundaryConditionsUpdater = formulation::updater::BoundaryConditionsUpdaterMock;

  auto single_group_test_ptr = dynamic_cast<GroupSolver*>(
      this->test_iterator_ptr_->group_solver_ptr());
  auto convergence_checker_test_ptr = dynamic_cast<ConvergenceChecker*>(
      this->test_iterator_ptr_->convergence_checker_ptr());
  auto moment_calculator_test_ptr = dynamic_cast<MomentCalculator*>(
      this->test_iterator_ptr_->moment_calculator_ptr());
  auto moment_map_test_ptr = dynamic_cast<MomentMapConvergenceChecker*>(
      this->test_iterator_ptr_->moment_map_convergence_checker_ptr());
  auto source_updater_test_ptr = dynamic_cast<SourceUpdater*>(
      this->test_iterator_ptr_->source_updater_ptr());
  auto boundary_conditions_test_ptr = dynamic_cast<BoundaryConditionsUpdater*>(
      this->test_iterator_ptr_->boundary_conditions_updater_ptr());

  EXPECT_NE(nullptr, single_group_test_ptr);
  EXPECT_NE(nullptr, convergence_checker_test_ptr);
  EXPECT_NE(nullptr, moment_calculator_test_ptr);
  EXPECT_NE(nullptr, moment_map_test_ptr);
  EXPECT_EQ(2, this->group_solution_ptr_.use_count());
  EXPECT_EQ(this->group_solution_ptr_.get(),
            this->test_iterator_ptr_->group_solution_ptr().get());
  EXPECT_NE(nullptr, source_updater_test_ptr);
  EXPECT_NE(nullptr, this->test_iterator_ptr_->reporter_ptr());
  EXPECT_NE(nullptr, boundary_conditions_test_ptr);
}

TYPED_TEST(IterationGroupSourceIterationTest, ConstructorThrows) {
  for (int i = 0; i < 5; ++i) {
    auto group_solver_ptr = (i == 0) ? nullptr :
        std::make_unique<solver::group::SingleGroupSolverMock>();
    auto convergence_checker_ptr = (i == 1) ? nullptr :
        std::make_unique<convergence::FinalCheckerMock<system::moments::MomentVector>>();
    auto moment_calculator_ptr = (i == 2) ? nullptr :
        std::make_unique<quadrature::calculators::SphericalHarmonicMomentsMock>();
    auto group_solution_ptr = (i == 3) ? nullptr :
        this->group_solution_ptr_;
    auto source_updater_ptr = (i == 4) ? nullptr :
        this->source_updater_ptr_;
    EXPECT_ANY_THROW({
      iteration::group::GroupSourceIteration<this->dim> test_iteration(
          std::move(group_solver_ptr),
          std::move(convergence_checker_ptr),
          std::move(moment_calculator_ptr),
          group_solution_ptr,
          source_updater_ptr
          );
    });
  }
}

TYPED_TEST(IterationGroupSourceIterationTest, ConstructorThrowNoBoundaryUpdater) {
  using BoundaryConditionsUpdater = formulation::updater::BoundaryConditionsUpdaterMock;

  auto group_solver_ptr = std::make_unique<solver::group::SingleGroupSolverMock>();
  auto convergence_checker_ptr = std::make_unique<convergence::FinalCheckerMock<system::moments::MomentVector>>();
  auto moment_calculator_ptr = std::make_unique<quadrature::calculators::SphericalHarmonicMomentsMock>();
  auto group_solution_ptr = this->group_solution_ptr_;
  auto source_updater_ptr = this->source_updater_ptr_;
  std::shared_ptr<BoundaryConditionsUpdater> boundary_condition_updater_ptr;

  EXPECT_ANY_THROW({
    iteration::group::GroupSourceIteration<this->dim> test_iteration(
        std::move(group_solver_ptr),
        std::move(convergence_checker_ptr),
        std::move(moment_calculator_ptr),
        group_solution_ptr,
        source_updater_ptr,
        boundary_condition_updater_ptr);
  });
}


template <typename DimensionWrapper>
class IterationGroupSourceSystemSolvingTest :
    public IterationGroupSourceIterationTest<DimensionWrapper> {
 public:

  using GroupSolution = typename IterationGroupSourceIterationTest<DimensionWrapper>::GroupSolution;

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
  std::array<dealii::PETScWrappers::MPI::Vector, total_groups> group_solutions_,
      group_rhs_;

  dealii::PETScWrappers::MPI::Vector expected_stored_solution_;


  int iterations = 0;
  const int max_iterations{1000};

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

    for (int angle = 0; angle < this->total_angles; ++angle) {
      this->energy_group_angular_solution_ptr_map_.insert(
          {system::SolutionIndex(group, angle),
           std::make_shared<dealii::Vector<double>>()});
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
  expected_stored_solution_.reinit(MPI_COMM_WORLD, 4, 4);
  for (int i = 0; i < 4; ++i)
    expected_stored_solution_[i] = 2.0;
  expected_stored_solution_.compress(dealii::VectorOperation::insert);
}

ACTION_P(Solve, test_class) {
  dealii::PETScWrappers::PreconditionNone no_conditioner(test_class->L_);
  test_class->solver_.solve(test_class->L_,
                            test_class->group_solutions_.at(arg0),
                            test_class->group_rhs_.at(arg0),
                            no_conditioner);
}

ACTION_P(Update, test_class) {
  auto& rhs = test_class->group_rhs_.at(arg1.get());
  auto& solution = test_class->group_solutions_.at(arg1.get());
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
  status.iteration_number = ++test_class->iterations;
  auto diff = arg0;
  diff.add(-1, arg1);
  status.delta = diff.l1_norm();
  if (status.delta.value() < 1e-7 || test_class->iterations >= test_class->max_iterations)
    status.is_complete = true;
  return status;
}

ACTION_P(ResetIterations, test_class) {
  test_class->iterations = 0;
}

TYPED_TEST_CASE(IterationGroupSourceSystemSolvingTest, bart::testing::AllDimensions);

TYPED_TEST(IterationGroupSourceSystemSolvingTest, UpdateThisAngularSolution) {
  this->test_iterator_ptr_->UpdateThisAngularSolutionMap(
      this->energy_group_angular_solution_ptr_map_);
  EXPECT_TRUE(this->test_iterator_ptr_->is_storing_angular_solution());
  EXPECT_EQ(this->test_iterator_ptr_->angular_solution_ptr_map().size(),
            this->total_groups*this->total_angles);
}

TYPED_TEST(IterationGroupSourceSystemSolvingTest, Iterate) {
  // This is the mock map to hold system current_moments
  using MockSolutionType = bart::system::solution::MPIGroupAngularSolutionMock;
  system::moments::MomentsMap current_moments, previous_moments;
  for (int group = 0; group < this->total_groups; ++group) {
    for (int l = 0; l <= this->max_harmonic_l; ++l) {
      for (int m = -l; m <= l; ++m) {
        system::moments::MomentIndex index{group, l, m};
        current_moments.emplace(index, 4);
        previous_moments.emplace(index, 4);
        EXPECT_CALL(*this->moments_obs_ptr_, BracketOp(index))
            .Times(AtLeast(1))
            .WillRepeatedly(ReturnRef(current_moments.at(index)));
        const auto& const_mock_current_moments = *this->moments_obs_ptr_;
        EXPECT_CALL(const_mock_current_moments, BracketOp(index))
            .Times(AtLeast(1))
            .WillRepeatedly(ReturnRef(current_moments.at(index)));
        EXPECT_CALL(*this->previous_moments_obs_ptr_, BracketOp(index))
            .Times(AtLeast(1))
            .WillRepeatedly(ReturnRef(previous_moments.at(index)));

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
      EXPECT_CALL(*this->source_updater_ptr_, UpdateScatteringSource(
          Ref(this->test_system),
          bart::system::EnergyGroup(group),
          quadrature::QuadraturePointIndex(angle)))
          .Times(AtLeast(1))
          .WillRepeatedly(Update(this));
      EXPECT_CALL(*this->boundary_conditions_updater_ptr_,
                  UpdateBoundaryConditions(
                      Ref(this->test_system),
                          bart::system::EnergyGroup(group),
                          quadrature::QuadraturePointIndex(angle)));
    }
  }

  auto mock_group_solution_ptr = dynamic_cast<MockSolutionType*>(this->group_solution_ptr_.get());

  convergence::Status moment_map_status;
  moment_map_status.is_complete = true;
  EXPECT_CALL(*this->moments_obs_ptr_, moments())
      .WillOnce(ReturnRef(current_moments));
  EXPECT_CALL(*this->moment_map_convergence_checker_obs_ptr_,
              CheckFinalConvergence(Ref(current_moments), _))
              .WillOnce(Return(moment_map_status));
  EXPECT_CALL(*this->moment_map_convergence_checker_obs_ptr_, Reset())
      .Times(AtLeast(1));

  EXPECT_CALL(*mock_group_solution_ptr, GetSolution(_))
      .Times(this->total_groups * this->total_angles)
      .WillRepeatedly(ReturnRef(this->expected_stored_solution_));

  EXPECT_CALL(*this->convergence_checker_obs_ptr_, Reset())
      .Times(AtLeast(1))
      .WillRepeatedly(ResetIterations(this));

  EXPECT_CALL(*this->convergence_checker_obs_ptr_, CheckFinalConvergence(_, _))
      .Times(AtLeast(1))
      .WillRepeatedly(ReturnConvergence(this));

  EXPECT_CALL(*this->reporter_ptr_, Report(A<const convergence::Status&>()))
      .Times(AtLeast(1));
  EXPECT_CALL(*this->reporter_ptr_, Report(A<const std::string&>()))
      .Times(AtLeast(1));

  this->test_system.total_groups = this->total_groups;
  this->test_system.total_angles = this->total_angles;
  
  EXPECT_CALL(*this->moments_obs_ptr_, max_harmonic_l())
      .WillRepeatedly(Return(this->max_harmonic_l));

  this->test_iterator_ptr_->UpdateThisAngularSolutionMap(
      this->energy_group_angular_solution_ptr_map_);

  this->test_iterator_ptr_->Iterate(this->test_system);
  EXPECT_LT(this->iterations, this->max_iterations);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(current_moments.at({0, 0, 0})[i],
                     this->true_scalar_flux_[i], 1e-6);
  }
  for (const auto& [index, solution_ptr] : this->energy_group_angular_solution_ptr_map_) {
    dealii::Vector<double> expected_solution;
    expected_solution = this->expected_stored_solution_;
    EXPECT_TRUE(test_helpers::CompareVector(expected_solution, *solution_ptr));
  }
}
} // namespace