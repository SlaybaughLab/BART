#ifndef BART_SRC_DOMAIN_DOMAIN_TYPES_H_
#define BART_SRC_DOMAIN_DOMAIN_TYPES_H_

#include <deal.II/dofs/dof_accessor.h>

namespace bart {

namespace domain {

template <int dim>
using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DOMAIN_TYPES_H_