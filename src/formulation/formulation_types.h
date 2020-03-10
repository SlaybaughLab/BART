#ifndef BART_SRC_FORMULATION_FORMULATION_TYPES_H_
#define BART_SRC_FORMULATION_FORMULATION_TYPES_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

namespace bart {

namespace formulation {

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

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_FORMULATION_TYPES_H_
