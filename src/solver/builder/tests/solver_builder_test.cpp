#include "solver/builder/solver_builder.hpp"

#include "solver/group/single_group_solver.h"
#include "solver/linear/gmres.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace solver = bart::solver;
namespace builder = solver::builder;
namespace test_helpers = bart::test_helpers;
using SolverName = builder::SolverName;

class SolverBuilderDefaultGMRESTest : public ::testing::Test {
 public:
  using ExpectedGroupSolver = solver::group::SingleGroupSolver;
  using ExpectedLinearSolver = solver::linear::GMRES;
  static constexpr int default_max_iterations{ 100 };
  static constexpr double default_convergence_tolerance { 1e-10 };
};

TEST_F(SolverBuilderDefaultGMRESTest, DefaultParameters) {
  auto solver_ptr = builder::SolverBuilder::BuildSolver(SolverName::kDefaultGMRESGroupSolver);
  ASSERT_NE(solver_ptr, nullptr);
  auto group_solver_ptr = dynamic_cast<ExpectedGroupSolver*>(solver_ptr.get());
  ASSERT_NE(group_solver_ptr, nullptr);
  auto linear_solver_ptr = dynamic_cast<ExpectedLinearSolver*>(group_solver_ptr->linear_solver_ptr());
  ASSERT_NE(linear_solver_ptr, nullptr);
  EXPECT_EQ(linear_solver_ptr->max_iterations(), default_max_iterations);
  EXPECT_EQ(linear_solver_ptr->convergence_tolerance(), default_convergence_tolerance);
}

TEST_F(SolverBuilderDefaultGMRESTest, SetParameters) {
  const int max_iterations { test_helpers::RandomInt(150, 200) };
  const double convergence_tolerance { test_helpers::RandomDouble(1e-10, 1e-6) };
  auto solver_ptr = builder::SolverBuilder::BuildSolver(SolverName::kDefaultGMRESGroupSolver,
                                                        max_iterations, convergence_tolerance);
  ASSERT_NE(solver_ptr, nullptr);
  auto group_solver_ptr = dynamic_cast<ExpectedGroupSolver*>(solver_ptr.get());
  ASSERT_NE(group_solver_ptr, nullptr);
  auto linear_solver_ptr = dynamic_cast<ExpectedLinearSolver*>(group_solver_ptr->linear_solver_ptr());
  ASSERT_NE(linear_solver_ptr, nullptr);
  EXPECT_EQ(linear_solver_ptr->max_iterations(), max_iterations);
  EXPECT_EQ(linear_solver_ptr->convergence_tolerance(), convergence_tolerance);
}

} // namespace
