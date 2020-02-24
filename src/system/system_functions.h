#ifndef BART_SRC_SYSTEM_SYSTEM_FUNCTIONS_H_
#define BART_SRC_SYSTEM_SYSTEM_FUNCTIONS_H_

#include "system/solution/mpi_group_angular_solution_i.h"
#include "domain/definition_i.h"

namespace bart {

namespace system {

/*! \brief Initializes all solutions and sets to a given value */
template <int dim>
void SetUpMPIAngularSolution(
    system::solution::MPIGroupAngularSolutionI &to_initialize,
    const domain::DefinitionI<dim> &domain_definition,
    const double value_to_set = 1.0);

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SYSTEM_FUNCTIONS_H_
