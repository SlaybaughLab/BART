#ifndef BART_SRC_FORMULATION_FORMULATION_TYPES_H_
#define BART_SRC_FORMULATION_FORMULATION_TYPES_H_

namespace bart {

namespace formulation {

template <int dim>
using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;

using Matrix = dealii::FullMatrix<double>;

enum class BoundaryType {
  kVacuum = 0,
  kReflective = 1
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_FORMULATION_TYPES_H_
