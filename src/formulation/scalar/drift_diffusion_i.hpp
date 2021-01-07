#ifndef BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include "domain/domain_types.h"
#include "system/system_types.h"

namespace bart::formulation::scalar {

template <int dim>
class DriftDiffusionI {
 public:
  using CellPtr = typename domain::CellPtr<dim>;
  using EnergyGroup = system::EnergyGroup;
  using Matrix = typename dealii::FullMatrix<double>;
  using Vector = typename dealii::Vector<double>;
  virtual ~DriftDiffusionI() = default;
  virtual auto FillCellDriftDiffusionTerm(Matrix& to_fill, const CellPtr&, const system::EnergyGroup,
                                          const Vector& group_scalar_flux,
                                          const Vector& integrated_angular_flux) const -> void = 0;
  virtual auto FillCellDriftDiffusionTerm(Matrix& to_fill, const CellPtr&, const system::EnergyGroup,
                                          const Vector& group_scalar_flux,
                                          const std::array<Vector, dim>& current) const -> void = 0;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_
