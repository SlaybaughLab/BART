#include "system/system_functions.h"

#include "system/terms/term.h"
#include "system/moments/spherical_harmonic.h"
#include "system/solution/mpi_group_angular_solution.h"

namespace bart {

namespace system {

template <int dim>
void SetUpMPIAngularSolution(
    system::solution::MPIGroupAngularSolutionI& to_initialize,
    const domain::DefinitionI<dim>& domain_definition,
    const double value_to_set) {
  auto& solution_map = to_initialize.solutions();
  AssertThrow(static_cast<int>(solution_map.size()) == to_initialize.total_angles(),
      dealii::ExcMessage("Error in SetUpMPIAngularSolution, MPIGroupAngularSolution solution map size does not match total_angles"))
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

void InitializeSystem(system::System &system_to_setup,
                 const int total_groups,
                 const int total_angles,
                 const bool is_eigenvalue_problem,
                 const bool is_rhs_boundary_term_variable) {
  using VariableLinearTerms = system::terms::VariableLinearTerms;

  std::string error_start{"Error: attempting to call Initialize System on a "
                          "system that appears to be already initialized: "};
  AssertThrow(system_to_setup.right_hand_side_ptr_ == nullptr,
              dealii::ExcMessage(error_start + "right hand side pointer is not null"))
  AssertThrow(system_to_setup.left_hand_side_ptr_ == nullptr,
              dealii::ExcMessage(error_start + "left hand side pointer is not null"))
  AssertThrow(system_to_setup.current_moments == nullptr,
              dealii::ExcMessage(error_start + "current moments pointer is not null"))
  AssertThrow(system_to_setup.previous_moments == nullptr,
              dealii::ExcMessage(error_start + "previous moments pointer is not null"))

  system_to_setup.total_groups = total_groups;
  system_to_setup.total_angles = total_angles;

  std::unordered_set<VariableLinearTerms> rhs_variable_terms{
    VariableLinearTerms::kScatteringSource};

  if (is_eigenvalue_problem) {
    system_to_setup.k_effective = 1.0;
    rhs_variable_terms.insert(VariableLinearTerms::kFissionSource);
  }

  if (is_rhs_boundary_term_variable) {
    rhs_variable_terms.insert(VariableLinearTerms::kReflectiveBoundaryCondition);
  }

  system_to_setup.right_hand_side_ptr_ = std::move(
      std::make_unique<system::terms::MPILinearTerm>(rhs_variable_terms));
  system_to_setup.left_hand_side_ptr_ = std::move(
      std::make_unique<system::terms::MPIBilinearTerm>());
  system_to_setup.current_moments = std::move(
      std::make_unique<system::moments::SphericalHarmonic>(total_groups, 0));
  system_to_setup.previous_moments = std::move(
      std::make_unique<system::moments::SphericalHarmonic>(total_groups, 0));
}

template <int dim>
void SetUpSystemTerms(system::System& system_to_setup,
                      const domain::DefinitionI<dim>& domain_definition) {
  const auto variable_terms =
      system_to_setup.right_hand_side_ptr_->GetVariableTerms();
  const int total_groups = system_to_setup.total_groups;
  const int total_angles = system_to_setup.total_angles;

  for (int group = 0; group < total_groups; ++group) {
    for (int angle = 0; angle < total_angles; ++angle) {
      system::Index index{group, angle};
      auto& lhs = system_to_setup.left_hand_side_ptr_;
      auto& rhs = system_to_setup.right_hand_side_ptr_;

      lhs->SetFixedTermPtr(index, domain_definition.MakeSystemMatrix());
      rhs->SetFixedTermPtr(index, domain_definition.MakeSystemVector());

      for (const auto variable_term : variable_terms) {
        rhs->SetVariableTermPtr(index, variable_term,
                                domain_definition.MakeSystemVector());
      }
    }
  }
}

void SetUpSystemMoments(system::System& system_to_setup,
                        const std::size_t solution_size) {

  auto initialize_moments = [=](system::moments::SphericalHarmonicI& to_initialize) {
    for (auto& moment : to_initialize) {
      moment.second.reinit(solution_size);
      moment.second = 1;
    }
  };

  initialize_moments(*system_to_setup.current_moments);
  initialize_moments(*system_to_setup.previous_moments);
}
void SetUpEnergyGroupToAngularSolutionPtrMap(
    solution::EnergyGroupToAngularSolutionPtrMap& to_setup,
    const int total_groups,
    const int total_angles) {
  using SolutionType = system::solution::MPIGroupAngularSolution;

  for (int group = 0; group < total_groups; ++group) {
    for (int angle = 0; angle < total_angles; ++angle) {
      to_setup.insert({system::SolutionIndex(group, angle),
                       std::make_shared<dealii::Vector<double>>()});
    };
  }
}

template void SetUpMPIAngularSolution<1>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<1>&, const double);
template void SetUpMPIAngularSolution<2>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<2>&, const double);
template void SetUpMPIAngularSolution<3>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<3>&, const double);

template void SetUpSystemTerms(system::System&, const domain::DefinitionI<1>&);
template void SetUpSystemTerms(system::System&, const domain::DefinitionI<2>&);
template void SetUpSystemTerms(system::System&, const domain::DefinitionI<3>&);

} // namespace system

} // namespace bart
