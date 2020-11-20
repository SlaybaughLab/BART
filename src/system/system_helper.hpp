#ifndef BART_SRC_SYSTEM_SYSTEM_HELPER_HPP_
#define BART_SRC_SYSTEM_SYSTEM_HELPER_HPP_

#include "domain/definition_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/system_helper_i.hpp"

namespace bart::system {

template <int dim>
class SystemHelper : public SystemHelperI<dim> {
 public:
  auto SetUpMPIAngularSolution(system::solution::MPIGroupAngularSolutionI &to_initialize,
                               const domain::DefinitionI<dim> &domain_definition,
                               const double value_to_set = 1.0) const -> void override;
};

} // namespace bart::system

#endif //BART_SRC_SYSTEM_SYSTEM_HELPER_HPP_
