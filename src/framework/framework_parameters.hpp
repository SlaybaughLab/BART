#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_

#include "problem/parameter_types.h"
#include "data/cross_sections.h"
#include "utility/named_type.h"

#include <memory>
#include <optional>
#include <set>

namespace bart::framework {

struct FrameworkParameters {
  using PolynomialDegree = utility::NamedType<int, struct PolynomialDegreeStruct>;

  int                         neutron_energy_groups{1};
  problem::EquationType       equation_type{problem::EquationType::kDiffusion};
  std::set<problem::Boundary> reflective_boundaries{};

  // Solver structure
  std::optional<problem::EigenSolverType> eigen_solver_type{std::nullopt};
  problem::InGroupSolverType              group_solver_type{problem::InGroupSolverType::kSourceIteration};

  // Solver domain
  problem::DiscretizationType discretization_type{problem::DiscretizationType::kContinuousFEM};
  problem::CellFiniteElementType cell_finite_element_type{problem::CellFiniteElementType::kGaussian};
  PolynomialDegree polynomial_degree{2};

  // Optional shared framework parts
  template <typename Part> using OptionalSharedPart = std::optional<std::shared_ptr<Part>>;
  OptionalSharedPart<data::CrossSections> cross_sections_;
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
