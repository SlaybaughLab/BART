#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_

#include "problem/parameter_types.h"

#include <set>

namespace bart::framework {

struct FrameworkParameters {
  int                         neutron_energy_groups{1};
  problem::EquationType       equation_type{problem::EquationType::kDiffusion};
  std::set<problem::Boundary> reflective_boundaries{};
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
