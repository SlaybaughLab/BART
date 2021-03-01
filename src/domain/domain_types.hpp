#ifndef BART_SRC_DOMAIN_DOMAIN_TYPES_HPP_
#define BART_SRC_DOMAIN_DOMAIN_TYPES_HPP_

#include <deal.II/dofs/dof_accessor.h>

#include "utility/named_type.h"

namespace bart::domain {

template <int dim>
using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;

using FaceIndex = bart::utility::NamedType<int, struct FaceIndexParameter>;

} // namespace bart::domain

#endif // BART_SRC_DOMAIN_DOMAIN_TYPES_HPP_