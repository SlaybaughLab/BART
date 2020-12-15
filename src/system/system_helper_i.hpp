#ifndef BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_
#define BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_

#include "domain/definition_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/solution/solution_types.h"
#include "system/system.hpp"


namespace bart::system {

template <int dim>
class SystemHelperI {
 public:
  virtual auto InitializeSystem(system::System& system_to_setup, const int total_groups,
                                const int total_angles, const bool is_eigenvalue_problem,
                                const bool is_rhs_boundary_term_variable) const -> void = 0;
  virtual auto SetUpEnergyGroupToAngularSolutionPtrMap(solution::EnergyGroupToAngularSolutionPtrMap& to_setup,
                                                       const int total_groups,
                                                       const int total_angles) const -> void = 0;
  virtual auto SetUpMPIAngularSolution(system::solution::MPIGroupAngularSolutionI &to_initialize,
                                       const domain::DefinitionI<dim> &domain_definition,
                                       const double value_to_set) const -> void = 0;
  virtual auto SetUpSystemTerms(system::System& system_to_setup,
                                const domain::DefinitionI<dim>& domain_definition) const -> void = 0;
  virtual auto SetUpSystemMoments(system::System& system_to_setup, const std::size_t solution_size) const -> void = 0;
};

} // namespace bart::system

#endif //BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_
