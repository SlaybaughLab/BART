#include "solver/builder/solver_builder.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace solver = bart::solver;
namespace builder = solver::builder;
using SolverName = builder::SolverName;


TEST(SolverBuilderTests, DefaultGroupSolverWithGMRES) {
  auto solver_ptr = builder::SolverBuilder::BuildSolver(SolverName::kDefaultGMRESGroupSolver);
  ASSERT_NE(solver_ptr, nullptr);
}

} // namespace
