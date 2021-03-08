#ifndef BART_SRC_FORMULATION_FORMULATION_TYPES_HPP_
#define BART_SRC_FORMULATION_FORMULATION_TYPES_HPP_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

namespace bart::formulation {

using FullMatrix = dealii::FullMatrix<double>;
using Vector = dealii::Vector<double>;

enum class BoundaryType {
  kVacuum = 0,
  kReflective = 1
};

// Formulation Implementations

enum class DiffusionFormulationImpl {
  kDefault = 0,
};

enum class SAAFFormulationImpl {
  kDefault = 0,
};

enum class StamperImpl {
  kDefault = 0,
};

} // bart::formulation

#endif //BART_SRC_FORMULATION_FORMULATION_TYPES_HPP_
