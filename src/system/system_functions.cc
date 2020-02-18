#include "system/system_functions.h"

namespace bart {

namespace system {

template <int dim>
void SetUpMPIAngularSolution(
    system::solution::MPIGroupAngularSolutionI& /*to_initialize*/,
    const domain::DefinitionI<dim>& /*domain_definition*/,
    const double /*value_to_set*/) {

}

template void SetUpMPIAngularSolution<1>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<1>&, const double);
template void SetUpMPIAngularSolution<2>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<2>&, const double);
template void SetUpMPIAngularSolution<3>(system::solution::MPIGroupAngularSolutionI&, const domain::DefinitionI<3>&, const double);


} // namespace system

} // namespace bart
