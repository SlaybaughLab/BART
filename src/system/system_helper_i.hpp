#ifndef BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_
#define BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_

namespace bart::system {

template <int dim>
class SystemHelperI {
 public:
  virtual auto SetUpMPIAngularSolution(system::solution::MPIGroupAngularSolutionI &to_initialize,
                                       const domain::DefinitionI<dim> &domain_definition,
                                       const double value_to_set) const -> void = 0;
};

} // namespace bart::system

#endif //BART_SRC_SYSTEM_SYSTEM_HELPER_I_HPP_
