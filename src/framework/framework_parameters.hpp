#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_

#include "problem/parameter_types.h"

#include <optional>
#include <set>

namespace bart::framework {

struct FrameworkParameters {
  int                         neutron_energy_groups{1};
  problem::EquationType       equation_type{problem::EquationType::kDiffusion};
  std::set<problem::Boundary> reflective_boundaries{};

  // Solver structure
  std::optional<problem::EigenSolverType> eigen_solver_type{std::nullopt};
  problem::InGroupSolverType              group_solver_type{problem::InGroupSolverType::kSourceIteration};

  // Solver domain
  problem::DiscretizationType discretization_type{problem::DiscretizationType::kContinuousFEM};
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
