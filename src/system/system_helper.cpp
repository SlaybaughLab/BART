#include "system/system_helper.hpp"

namespace bart::system {

template<int dim>
auto SystemHelper<dim>::SetUpMPIAngularSolution(solution::MPIGroupAngularSolutionI &to_initialize,
                                                const domain::DefinitionI<dim> &domain_definition,
                                                const double value_to_set) const -> void {
  auto& solution_map = to_initialize.solutions();
  AssertThrow(static_cast<int>(solution_map.size()) == to_initialize.total_angles(),
              dealii::ExcMessage("Error in SetUpMPIAngularSolution, MPIGroupAngularSolution solution map size does not "
                                 "match total_angles"))
  const auto locally_owned_dofs = domain_definition.locally_owned_dofs();

  for (auto& solution_pair : solution_map) {
    auto& solution = solution_pair.second;
    solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    auto local_elements = solution.locally_owned_elements();
    for (auto index : local_elements) {
      solution[index] = value_to_set;
    }
    solution.compress(dealii::VectorOperation::insert);
  }
}

template class SystemHelper<1>;
template class SystemHelper<2>;
template class SystemHelper<3>;
} // namespace bart::system


