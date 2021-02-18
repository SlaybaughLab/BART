#ifndef BART_SRC_FORMULATION_COMMON_RHS_CONSTANT_I_HPP_
#define BART_SRC_FORMULATION_COMMON_RHS_CONSTANT_I_HPP_

#include "domain/domain_types.h"

#include <deal.II/lac/vector.h>

namespace bart::formulation::common {

template <int dim>
class RHSConstantI {
 public:
  using CellPtr = domain::CellPtr<dim>;
  using Vector = dealii::Vector<double>;
  virtual ~RHSConstantI() = default;
  virtual auto FillCellConstantTerm(Vector& to_fill, const CellPtr&, const Vector& constant_vector) const -> void = 0;
};

} // namespace bart::formulation::common

#endif //BART_SRC_FORMULATION_COMMON_RHS_CONSTANT_I_HPP_
